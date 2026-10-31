/* Included only in the isolated RS source copy, after rp_setup definitions. */
static int panel_enabled=0,panel_h_step=0;
static uint64_t panel_revision=0,panel_h_hash=0,panel_s_hash=0,panel_geometry=0;
static double ****panel_h_field=NULL,****panel_s_field=NULL;
static uint64_t panel_hash_bytes(uint64_t h,const void*data,size_t bytes){const unsigned char*x=data;for(size_t i=0;i<bytes;i++){h^=x[i];h*=UINT64_C(1099511628211);}return h;}
static uint64_t panel_geometry_now(void){uint64_t h=UINT64_C(1469598103934665603);for(int a=1;a<=3;a++)h=panel_hash_bytes(h,tv[a]+1,sizeof(double)*3);for(int a=1;a<=atomnum;a++)h=panel_hash_bytes(h,Gxyz[a]+1,sizeof(double)*3);return h?h:1;}
static uint64_t panel_field_hash(double****field,double*packed){
 uint64_t hash=UINT64_C(1469598103934665603);size_t at=0;
 for(int q=0;q<rp_nr;q++){RSBlock*b=rs_block+rp_id[q];int mc=F_G2M[b->a];for(int i=0;i<b->na;i++)for(int j=0;j<b->nb;j++){double v=field[mc][b->h][i][j];if(!isfinite(v))rs_die("nonfinite representative projected field seal");if(packed)packed[at]=v;at++;hash=panel_hash_bytes(hash,&v,sizeof(v));}}
 if(at!=rp_entries)rs_die("representative projected seal order mismatch");return hash;
}
int CpuRSym_PanelEnable(int enabled){
 int lo,hi;MPI_Allreduce(&enabled,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(&enabled,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);
 if(lo!=hi||(enabled!=0&&enabled!=1)||panel_revision)rs_die("invalid panel enable agreement/lifecycle");panel_enabled=enabled;return panel_enabled;
}
static void cpu_panel_record_projected(int step,double****h,double****s){
 if(!panel_enabled)return;if(CpuRSym_Enabled()!=2||step<1||step<=panel_h_step)rs_die("panel seal needs fresh per-SCF Reynolds projection");rp_setup();if(rs_ng!=6)rs_die("representative sector assembly currently requires D3");
 uint64_t geometry=panel_geometry_now();if(panel_geometry&&geometry!=panel_geometry)rs_die("panel seal fixed geometry changed");panel_geometry=geometry;
 if(!panel_s_field){panel_s_field=s;panel_s_hash=panel_field_hash(s,NULL);}else if(panel_s_field!=s||panel_s_hash!=panel_field_hash(s,NULL))rs_die("fixed projected overlap changed after seal");
 panel_h_field=h;panel_h_hash=panel_field_hash(h,NULL);panel_h_step=step;panel_revision++;
}
int CpuRSym_PanelAOMetadata(int*nat,int*species,int*ng,const int**types,const int**norb,const int*const**permutations,const double*const**rotations){
 if(!panel_enabled||!rs_ready||!panel_revision||rs_ng!=6)return 0;rs_geometry_guard();*nat=rs_nat;*species=rs_ns;*ng=rs_ng;*types=rs_type;*norb=rs_norb;*permutations=(const int*const*)rs_perm;*rotations=(const double*const*)rs_d;return 1;
}
int CpuRSym_PanelFields(double****h,double****s){return panel_enabled&&panel_revision&&h==panel_h_field&&s==panel_s_field;}
RpaBlock *CpuRSym_PanelBlocks(int*count){
 if(!panel_enabled||!panel_revision)rs_die("representative panel metadata requested before projection seal");rp_setup();int*ao=rs_alloc(rs_nat+1,sizeof(int));for(int a=0;a<rs_nat;a++)ao[a+1]=ao[a]+rs_norb[rs_type[a]];
 RpaBlock*out=rs_alloc(rp_nr,sizeof(*out));for(int q=0;q<rp_nr;q++){RSBlock*b=rs_block+rp_id[q];out[q]=(RpaBlock){ao[b->a-1],b->na,ao[b->b-1],b->nb,rp_mult[q],(size_t)rp_offsets[q]};}free(ao);*count=rp_nr;return out;
}
RpaContext CpuRSym_PanelContext(uint64_t basis,uint64_t distribution){return(RpaContext){panel_geometry,basis,distribution,rs_map_fnv};}
int CpuRSym_PanelPack(int kind,int step,RpaContext context,double*values,size_t entries,RpaStamp*stamp){
 if(!panel_enabled||!panel_revision||entries!=rp_entries||(!values&&entries)||!stamp||context.geometry!=panel_geometry||context.geometry!=panel_geometry_now()||context.projector_map!=rs_map_fnv||!context.basis||!context.distribution)return 0;
 double****field=kind==1?panel_h_field:kind==2?panel_s_field:NULL;if(!field||(kind==1&&step!=panel_h_step)||(kind==2&&step!=0))return 0;
 uint64_t seal=kind==1?panel_h_hash:panel_s_hash,actual=panel_field_hash(field,values);if(actual!=seal)return 0;
 *stamp=(RpaStamp){context,kind==1?panel_revision:1,actual,kind,step};return 1;
}
