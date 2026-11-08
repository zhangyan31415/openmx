"""Replot the saved cutoff comparison; requires NumPy and Matplotlib.

Run `python plot.py` in any directory. Input CSVs are resolved beside this file.
No interpolation, new SCF calculation or per-band alignment is performed.
"""
from pathlib import Path
import csv
import json
import hashlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parent
BLUE = 'dodgerblue'
RED = '#DE6E66'


def style():
    # Use an available local serif font; no bundled or required proprietary fonts.
    fonts = Path.home() / '.local/share/fonts'
    for name in ('times.ttf', 'timesbd.ttf', 'timesi.ttf', 'timesbi.ttf'):
        path = fonts / name
        if path.is_file():
            font_manager.fontManager.addfont(str(path))
    available = {f.name for f in font_manager.fontManager.ttflist}
    serif = 'Times New Roman' if 'Times New Roman' in available else 'DejaVu Serif'
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': [serif],
        'mathtext.fontset': 'stix', 'font.size': 11,
        'axes.linewidth': .7, 'axes.unicode_minus': False,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'xtick.major.size': 3, 'ytick.major.size': 3,
        'xtick.major.width': .6, 'ytick.major.width': .6,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white',
    })


def load(material):
    metadata = json.loads((ROOT / f'{material}.json').read_text())
    path = ROOT / f'{material}.csv'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == metadata['csv_sha256']
    rows = list(csv.DictReader(path.open()))
    old = np.empty((13, 8)); new = np.empty_like(old)
    k = np.empty((13, 2)); filled = np.zeros_like(old, dtype=bool)
    for row in rows:
        i, j = int(row['path_index']), int(row['band_from_top']) + 7
        assert not filled[i, j]
        old[i, j], new[i, j] = float(row['old_eV']), float(row['new_eV'])
        k[i] = float(row['kx']), float(row['ky'])
        filled[i, j] = True
    assert filled.all() and np.isfinite(old).all() and np.isfinite(new).all()
    shift = (new-old).mean()
    delta = 1000*(new-old-shift)
    metrics = metadata['metrics']
    np.testing.assert_allclose([np.sqrt((delta**2).mean()), np.abs(delta).max()],
        [metrics['aligned_RMS_difference_meV'], metrics['aligned_max_difference_meV']], atol=1e-7)
    # Hexagonal direct vectors enclose 120 degrees; reciprocal vectors enclose 60.
    dk = np.diff(k, axis=0)
    ds = np.sqrt(dk[:, 0]**2 + dk[:, 1]**2 + dk[:, 0]*dk[:, 1])
    x = np.r_[0, np.cumsum(ds)]; x /= x[-1]
    zero = old.max()
    return dict(name=material, x=x, old=1000*(old-zero),
        new=1000*(new-shift-zero), delta=delta, metrics=metrics)


def draw(records, name):
    ncol = len(records)
    fig, axes = plt.subplots(4, ncol, squeeze=False,
        figsize=(4.2*ncol+.45, 10.0),
        gridspec_kw={'height_ratios': [1, 1, 1, .65]})
    fig.subplots_adjust(left=.09 if ncol > 1 else .17, right=.975,
        bottom=.075, top=.865, wspace=.29, hspace=.38)
    fig.suptitle('Integration cutoff: valence-band comparison', fontsize=15, y=.985)
    fig.text(.5, .951, r'$2.13393^\circ$  ·  4326 atoms per material  ·  NSOC',
        ha='center', fontsize=11)
    fig.legend([Line2D([], [], color=BLUE, lw=1),
                Line2D([], [], color=RED, lw=1, ls=(0, (4, 2)))],
        [r'150 Ry / $10^{-7}$', r'80 Ry / $10^{-6}$ (aligned)'],
        loc='upper center', bbox_to_anchor=(.5, .937), ncol=2,
        frameon=False, handlelength=2.5, columnspacing=2.1, fontsize=10.5)
    floor = 10*np.floor(min(min(r['old'].min(), r['new'].min()) for r in records)/10)
    for col, rec in enumerate(records):
        x = rec['x']; ticks = x[[0, 4, 8, 12]]
        material = r'MoTe$_2$' if rec['name'] == 'MoTe2' else r'WSe$_2$'
        stages = ['Overlay', '150 Ry', '80 Ry (aligned)', 'Aligned difference']
        for row, ax in enumerate(axes[:, col]):
            for xpos in ticks[1:-1]: ax.axvline(xpos, color='.78', lw=.65, zorder=0)
            ax.axhline(0, color='.83', lw=.6, zorder=0)
            ax.set_xlim(0, 1)
            ax.set_xticks(ticks, [r'$\Gamma$', 'M', 'K', r'$\Gamma$'])
            ax.tick_params(labelsize=10, pad=3)
            label = f'({chr(97+row*ncol+col)}) '
            ax.set_title(label + (material+' · ' if row == 0 else '') + stages[row],
                loc='left', fontsize=11.5, pad=7)
            if row < 3:
                ax.set_ylim(floor-1, 2)
                ax.set_yticks(np.arange(0, floor-1, -20))
                if row in (0, 1): ax.plot(x, rec['old'], color=BLUE, lw=.85)
                if row in (0, 2): ax.plot(x, rec['new'], color=RED, lw=.9,
                    ls=(0, (4, 2)) if row == 0 else '-')
                if col == 0: ax.set_ylabel('Energy (meV)', labelpad=7)
            else:
                ax.plot(x, rec['delta'], color=RED, lw=.65, alpha=.85)
                ax.set_ylim(-.5, .5); ax.set_yticks([-.4, 0, .4])
                if col == 0: ax.set_ylabel(r'$\Delta E$ (meV)', labelpad=7)
            if row == 0:
                m = rec['metrics']
                ax.text(.97, .06, f"RMS {m['aligned_RMS_difference_meV']:.3f} meV\n"
                    f"Max {m['aligned_max_difference_meV']:.3f} meV",
                    transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=2))
    fig.text(.5, .025, r'TAPW $K_1$, Q-shell 8 · 13 k points · highest 8 valence bands',
        ha='center', fontsize=10, color='.25')
    for ext in ('png', 'pdf', 'svg'):
        output = ROOT/f'{name}.{ext}'
        fig.savefig(output, dpi=220, bbox_inches='tight', pad_inches=.07)
        if ext == 'svg':
            output.write_text('\n'.join(line.rstrip() for line in output.read_text().splitlines())+'\n')
    plt.close(fig)


def main():
    style(); records = [load(m) for m in ('MoTe2', 'WSe2')]
    draw(records, 'cutoff_comparison')
    for r in records: draw([r], r['name']+'_aligned')
    print('Figures regenerated; CSV hashes and aligned metrics verified.')


if __name__ == '__main__':
    main()
