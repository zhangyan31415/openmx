#!/usr/bin/env python3
"""Apply verified NSOC deltas onto SOC sources with exact, unique context.

No fuzzy merge: unknown edits fail before producing output. Callers retain the
stock source, archived NSOC source and SOC source hashes in their build manifest.
"""
import argparse
import hashlib
import json
import difflib
import re
from pathlib import Path


def merge_text(base: str, nsoc: str, soc: str) -> str:
    """Three-way line mapping; repeated native blocks cannot move an edit.

    Coordinates come from the original stock-to-SOC alignment, then edits are
    applied backwards. Context lookups in already-mutated text are forbidden.
    """
    if soc == base:
        return nsoc
    if base and soc.count(base)>1:
        raise ValueError('ambiguous duplicate stock context')
    original=base.splitlines(keepends=True)
    scalar=nsoc.splitlines(keepends=True)
    spin=soc.splitlines(keepends=True)
    mapping={}
    for block in difflib.SequenceMatcher(None,original,spin,autojunk=False).get_matching_blocks():
        for offset in range(block.size):
            mapping[block.a+offset]=block.b+offset
    edits=[]
    for tag,a,b,c,d in difflib.SequenceMatcher(None,original,scalar,autojunk=False).get_opcodes():
        if tag=='equal':continue
        new=scalar[c:d]
        previous=max((i for i in mapping if i<a),default=-1)
        following=min((i for i in mapping if i>=b),default=len(original))
        low=mapping[previous]+1 if previous>=0 else 0
        high=mapping[following] if following<len(original) else len(spin)
        if a==b:
            # A SOC insertion at this boundary is retained. A shared insertion
            # already present is recognized only within this mapped boundary.
            existing=''.join(spin[low:high])
            if ''.join(new) and ''.join(new) in existing:continue
            if previous!=a-1 or following!=a:
                raise ValueError(f'NSOC insertion overlaps edited stock lines {a+1}:{b}')
            edits.append((low,low,new))
        elif all(i in mapping for i in range(a,b)) and all(mapping[i]==mapping[a]+i-a for i in range(a,b)):
            edits.append((mapping[a],mapping[b-1]+1,new))
        elif spin[low:high]==new:
            continue
        else:
            raise ValueError(f'NSOC {tag} conflicts with SOC edit at stock lines {a+1}:{b}')
    for low,high,new in sorted(edits,key=lambda item:item[0],reverse=True):
        spin[low:high]=new
    return ''.join(spin)


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'expected unique native context: {old[:100]!r}')
    return text.replace(old, new, 1)


def merge_truncation(nsoc: str, soc: str) -> str:
    """Keep SOC allocation guards and add scalar-only history/ghost behavior.

    The two archived files share grid/field/VNA helper ABI. Their validity
    predicates differ intentionally; apply each within its original spin path.
    """
    imports = ('#include "cpu_nsoc_imag.h"\n'
               '#include "cpu_nsoc_symhistory_guard.h"\n'
               'extern void CpuNSocSymHistoryDestroy(void);\n'
               '#include "cpu_nsoc_history_capacity.h"\n'
               '#include "cpu_nsoc_ghost.h"\n')
    out = replace_once(soc, '#include "mpi.h"\n', '#include "mpi.h"\n' + imports)
    configure = ('  if(CpuGridFields_Configure(mpi_comm_level1,SpinP_switch==3?4:2)!=0){\n'
                 '    fprintf(stderr,"CPU_GRID_FIELDS_LAYOUT_REJECT inconsistent scalar/spinor layout\\n");\n'
                 '    MPI_Abort(mpi_comm_level1,97);exit(97);\n'
                 '  }\n')
    out = replace_once(out, '  MPI_Comm_rank(mpi_comm_level1,&myid);\n  {\n    int valid=',
                       '  MPI_Comm_rank(mpi_comm_level1,&myid);\n' + configure + '  {\n    int valid=')
    scalar_valid = re.search(r'    int valid=(.*?);', nsoc, re.S).group(1)
    original_valid = re.search(r'    int valid=(.*?);', soc, re.S).group(0)
    soc_valid = original_valid.split('=', 1)[1][:-1]
    out = replace_once(out, original_valid,
        '    int valid=(SpinP_switch==0 ? ('+scalar_valid+') : ('+soc_valid+'));')
    old = '  int cpu_local_hs=CpuLocalHSPlan(getenv("OPENMX_CPU_LOCAL_HS"),&local_hs_input);'
    new = ('  int cpu_local_hs=0,cpu_local_olp=0;\n'
           '  if(SpinP_switch==0){\n'
           '    if(CpuNSOCGhostAgree(mpi_comm_level1,('+scalar_valid+'),&cpu_local_hs,&cpu_local_olp)<0){\n'
           '      fprintf(stderr,"CPU_NSOC_GHOST_REJECT scalar owner-only flags require validated NSOC path\\n");\n'
           '      MPI_Abort(mpi_comm_level1,96);exit(96);\n'
           '    }\n'
           '  }else cpu_local_hs=CpuLocalHSPlan(getenv("OPENMX_CPU_LOCAL_HS"),&local_hs_input);')
    out = replace_once(out, old, new)
    old = '  int cpu_local_olp=CpuLocalOLPPlan(getenv("OPENMX_CPU_LOCAL_OLP"),getenv("OPENMX_CPU_LOCAL_HS"),&local_hs_input);'
    out = replace_once(out, old, '  if(SpinP_switch!=0)cpu_local_olp=CpuLocalOLPPlan(getenv("OPENMX_CPU_LOCAL_OLP"),getenv("OPENMX_CPU_LOCAL_HS"),&local_hs_input);')
    # Apply all history/imaginary-only deltas from the NSOC source. Filtering
    # by added API avoids touching shared VNA/SOC allocation changes.
    stock = Path(__file__).resolve().parents[1]/'vendor/openmx4.0/source/truncation.c'
    base_lines = stock.read_text().splitlines(keepends=True)
    nsoc_lines = nsoc.splitlines(keepends=True)
    for tag,a,b,c,d in difflib.SequenceMatcher(None,base_lines,nsoc_lines,autojunk=False).get_opcodes():
        if tag=='equal':
            continue
        old,new=''.join(base_lines[a:b]),''.join(nsoc_lines[c:d])
        if not any(api in new for api in ('CpuNSocHamiltonianHistorySlots()',
                'CpuNSocSymHistoryEnabled()', 'CpuNSOCImagSlots()', 'CpuNSOCImagReport(')):
            continue
        # Top-level include and guard blocks are handled above.
        if '#include' in new or 'int valid=' in new:
            continue
        applied=False
        for width in (3,6,12,1):
            before=''.join(base_lines[max(0,a-width):a])
            after=''.join(base_lines[b:b+width])
            context=before+old+after
            if out.count(context)==1:
                out=out.replace(context,before+new+after,1)
                applied=True
                break
        if not applied:
            raise ValueError(f'truncation history conflict stock lines {a+1}:{b}')
    return out


def make_unified_selection(repo: Path) -> dict:
    """Generate a single build selection with scalar and spinor paths retained."""
    soc = json.loads((repo/'cpu/overrides/soc/SOURCES.json').read_text())
    nsoc = json.loads((repo/'cpu/overrides/nsoc/SOURCES.json').read_text())
    by_name = {u['object'].split('.')[0]: u for u in soc['units']}
    replacements, extra, provenance = {}, [], []
    outputs = {}
    for unit in nsoc['units']:
        key = unit['object'].split('.')[0]
        source = repo/unit['source']
        if key not in by_name:
            extra.append(unit)
            continue
        spin = by_name[key]
        old = repo/spin['source']
        if source.read_bytes() == old.read_bytes():
            continue
        if key == 'Force':
            # Batch HNL/one-shot force extends the same shared-field native
            # force path; do not erase SOC optimizations with its predecessor.
            provenance.append({'object': spin['object'], 'policy': 'retain SOC batch force superset of native NSOC field path'})
            continue
        if key in ('Set_Hamiltonian','Set_Aden_Grid','Set_Orbitals_Grid','Set_dOrbitals_Grid'):
            # Flag-off path is the same native field assembly; opt-in BLAS3
            # is guarded for scalar NSOC and leaves spinor assembly untouched.
            replacements[spin['object']] = unit['source']
            continue
        if key == 'truncation':
            merged = merge_truncation(source.read_text(), old.read_text())
        else:
            stock = repo/'cpu/vendor/openmx4.0/source'/f'{key}.c'
            if not stock.exists():
                raise ValueError(f'no stock baseline for shared {key}')
            merged = merge_text(stock.read_text(), source.read_text(), old.read_text())
        if key == 'DFT':
            # LOCAL_HS is shared by scalar/full-SCF and spinor/one-shot
            # allocation plans. Only the latter has a one-iteration bound.
            merged = replace_once(merged,
                'if((CpuOwnerSparseRequested()||(local_hs&&!strcmp(local_hs,"1")))&&DFTSCF_loop!=1){',
                'if(SpinP_switch==3 && SO_switch==1 && (CpuOwnerSparseRequested()||(local_hs&&!strcmp(local_hs,"1")))&&DFTSCF_loop!=1){')
        # Later verified generic-memory header admits every positive input
        # history while retaining the original dc095 archive byte-for-byte.
        merged = merged.replace('#include \"cpu_nsoc_history_capacity.h\"', '#include \"cpu_nsoc_history_capacity.unified.h\"')
        destination = repo/'cpu/overrides/nsoc/source'/f'{key}.unified.c'
        outputs[destination] = merged
        replacements[spin['object']] = str(destination.relative_to(repo))
        provenance.append({'object': spin['object'], 'nsoc_source': unit['source'],
            'soc_source': spin['source'], 'sha256': hashlib.sha256(merged.encode()).hexdigest()})
    replacements['Set_XC_Grid.o'] = 'cpu/overrides/nsoc/source/Set_XC_Grid.workspace.c'
    replacements['cpu_grid_fields_shared.o'] = 'cpu/overrides/nsoc/source/cpu_grid_fields_shared.unified.c'
    # Source equality is insufficient for helpers whose private header used
    # different scalar/spinor layout constants. Every compiled consumer must
    # include the same explicit unified header, regardless of source family.
    for unit in soc['units']:
        selected = repo/replacements.get(unit['object'],unit['source'])
        text = outputs.get(selected)
        if text is None:
            text = selected.read_text()
        if '#include "cpu_grid_fields_shared.h"' not in text:
            continue
        text = text.replace('#include "cpu_grid_fields_shared.h"', '#include "cpu_grid_fields_shared.unified.h"')
        key = unit['object'].split('.')[0]
        destination = repo/'cpu/overrides/nsoc/source'/f'{key}.unified.c'
        outputs[destination] = text
        replacements[unit['object']] = str(destination.relative_to(repo))
    # All conflicts must be resolved before any generated file is overwritten.
    for destination, text in outputs.items():
        destination.write_text(text)
    return {'schema': 'openmx-unified-cpu-selection-v1', 'replacements': replacements,
            'extra_units': extra, 'include_dirs': nsoc['include_dirs'],
            'capabilities': {'nsoc_symmetry': True, 'nsoc_positive_history': True, 'nsoc_xc_workspace': True},
            'provenance': provenance}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stock', type=Path)
    parser.add_argument('--nsoc', type=Path)
    parser.add_argument('--soc', type=Path)
    parser.add_argument('--unified-repo', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.unified_repo:
        text = json.dumps(make_unified_selection(args.unified_repo.resolve()), indent=2) + '\n'
    else:
        if not all((args.stock, args.nsoc, args.soc)):
            parser.error('--stock, --nsoc and --soc are required without --unified-repo')
        text = merge_text(args.stock.read_text(), args.nsoc.read_text(), args.soc.read_text())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text)

if __name__ == '__main__':
    main()
