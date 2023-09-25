'''
Setup flopy config files 

One-dimensional column with equivalent dispersion binary reactive system
'''

# Dependencies
import os
import sys
import numpy as np
import flopy
import matplotlib.pyplot as plt
import argparse
import xmipy
# self
import config
import scipy.special as spc

##################################################################################################
# Arguments
parser = argparse.ArgumentParser( description='Setup simulations for experiment.' )
parser.add_argument( '--write'    , action='store_true', help='write simulation files.' )
parser.add_argument( '--run'      , action='store_true', help='executes simulation, requires mf6lib.' )
parser.add_argument( '--force'    , action='store_true', help='skip the prompt after verification of time step.' )
args   = parser.parse_args()


#########################
# General configuration #
#########################
stressperiods = [
        {
            'id'           : 0,
            'length'       : 3*10*24.0, # 30 days in hours 
            'n_time_steps' : 3*21599,   # 
            'ts_multiplier': 1,
            'steady_state' : False,
            'printfreqts'  : 1, 
        },
    ]
experiment_folder = config.FOLDERS['sims']
base_model_name   = 'mf6model'
SAVEALL           = True # needed if the integrated reaction product aims to be calculated

component = 0.5
dpdsdiff  = 0.1 # dp/ds
sim_name  = 'NLCOMP1D.dpds.'+str(dpdsdiff)+'.u'+str(component)
davg      = 1e-9*1e4*3600   # avg diffusion, the same for all sims
dmscale   = 1e-8*1e4*3600   # scale for time step selection 
dmeffsec  = 2*davg/( 1 + dpdsdiff )
dmeffprim = dpdsdiff*dmeffsec
sim_ws    = os.path.join( experiment_folder, sim_name )

####################
# Setup flow model #
####################

# sim 
sim = flopy.mf6.MFSimulation(
        sim_name=sim_name,
        sim_ws=sim_ws,
        exe_name='mf6'
    )

# tdis
perioddata = []
for sp in stressperiods:
    perioddata.append( [sp['length'], sp['n_time_steps'], sp['ts_multiplier']] )
tdis = flopy.mf6.ModflowTdis(
    sim,
    nper=len(stressperiods),
    perioddata=perioddata,
    time_units='hours',
)

# ims 
hclose     = 1e-7
ninner     = 100
nouter     = 50
complexity = 'COMPLEX'
accel      = 'BICGSTAB'
ims = flopy.mf6.ModflowIms(
        sim,
        pname         = 'ims',
        print_option  = 'SUMMARY',
        complexity    = complexity,
        outer_dvclose = hclose, 
        outer_maximum = nouter, 
        inner_maximum = ninner, 
        inner_dvclose = hclose, 
        linear_acceleration = accel,
    )

# gwf
gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=base_model_name,
        save_flows=True
    )

# dis
nlay = 1
nrow = 1
ncol = 500
delc = 100.0/ncol
delz = delc
delr = delc
dis  = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay, nrow=nrow, ncol=ncol,
        delr=delr, delc=delc, top=delz, botm=0,
        length_units='CENTIMETERS',
    )

# npf
hk  = 0.01
npf = flopy.mf6.ModflowGwfnpf(
    gwf,
    save_specific_discharge=True,
    icelltype=0,
    k=hk,
)

# sto
npf = flopy.mf6.ModflowGwfsto(gwf, ss=0,sy=0, steady_state=True)

# chd (no flow) 
hdin   = delc
hdout  = hdin
chdspd = []
for ir in range(nrow):
    chdspd.append( [ (0, ir,      0), hdin  ] )
    chdspd.append( [ (0, ir, ncol-1), hdout ] )
chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        maxbound=len(chdspd),
    )

# ic
flopy.mf6.ModflowGwfic(
    gwf,
    strt=hdin,
)

# oc
head_file   = base_model_name + '.hds'
budget_file = base_model_name + '.bud'   
if SAVEALL:
    saverecord=[ ('HEAD', 'ALL') ]
else:
    saverecord=[ ('HEAD', 'FIRST'), ('HEAD','LAST') ]
flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord  =head_file,
    budget_filerecord=budget_file,
    saverecord=saverecord
)

#########################
# Setup transport model #
#########################
if SAVEALL:
    gwtname = "gwtall-" + base_model_name
else:
    gwtname = "gwt-" + base_model_name

# gwt
gwt = flopy.mf6.MFModel(
    sim,
    model_type="gwt6",
    modelname=gwtname,
    model_nam_file="{}.nam".format(gwtname),
)
gwt.name_file.save_flows = True

#imsgwt
imsgwt = flopy.mf6.ModflowIms(
    sim,
    print_option        = "SUMMARY",
    outer_dvclose       = hclose,
    outer_maximum       = nouter,
    under_relaxation    = "NONE",
    inner_maximum       = ninner,
    inner_dvclose       = hclose,
    linear_acceleration = accel, 
    complexity          = complexity,
    filename="{}.ims".format(gwtname),
)
sim.register_ims_package(imsgwt, [gwt.name])

# disgwt
disgwt = flopy.mf6.ModflowGwfdis(
        gwt,
        nlay=nlay, nrow=nrow, ncol=ncol,
        delr=delr, delc=delc, top=delz, botm=0,
    )

# adv
scheme = "TVD"
flopy.mf6.ModflowGwtadv(
    gwt,
    scheme=scheme,
    filename="{}.adv".format(gwtname)
)

#ic (smoothed)
def analytical_erf( coords, diffusion, time, uzero, uend=0):
    erf_function = spc.erf( (coords)/np.sqrt(4*diffusion*time) )
    return 0.5*( uzero-uend)*(1-erf_function) + uend
xv  = np.arange(0,ncol,1)*delc
erf = analytical_erf( xv - 50, davg, 2, component, -1*component )
sconc = np.zeros((1,nrow,ncol))
sconc[0,0,:] = erf
##ic (sharp)
#sconc[0,:,0:int(0.5*ncol)] = component
#sconc[0,:,int(0.5*ncol): ] = -1*component
flopy.mf6.ModflowGwtic(
    gwt,
    strt=sconc,
    filename="{}.ic".format(gwtname)
)

# dsp
porosity    = 1
deltax      = delc 
tscale      = deltax**2/(dmscale)
ct          = 0.1
trec        = ct*tscale

print( ' # Dispersion analysis:  ' )
print( '   - diffp               : ', dmeffprim/porosity                              , ' cm2/hr' )
print( '   - diffs               : ', dmeffsec /porosity                              , ' cm2/hr' )
print( '   - diffp               : ', dmeffprim/porosity/1e4/3600                     , ' m2/s'   )
print( '   - diffs               : ', dmeffsec /porosity/1e4/3600                     , ' m2/s'   )
print( '   - diffavg             : ', np.mean([dmeffprim,dmeffsec])/porosity          , ' cm2/hr' )
print( '   - diffavg             : ', np.mean([dmeffprim,dmeffsec])/porosity/1e4/3600 , ' m2/s'   )
print( '   - dp/ds               : ', dpdsdiff                                                    )
print( '   - deltat (given)      : ', stressperiods[0]['length']/stressperiods[0]['n_time_steps'] )
print( '   - deltat (recommended): ', trec                                                        )
print( '   - nsteps (given)      : ', stressperiods[0]['n_time_steps']                            ) 
print( '   - nsteps (recommended): ', stressperiods[0]['length']/trec                             ) 
print( '   - effective CT        : ', (stressperiods[0]['length']/stressperiods[0]['n_time_steps'])/trec*ct  ) 
print( 'Continue ... ? (press c to continue, or q to exit)' )
if not args.force:
    import pdb
    pdb.set_trace()

nlfunc     = sconc/np.sqrt(sconc**2+1)
equivdiffc = 0.5*( ( dmeffprim  - dmeffsec  )*nlfunc + dmeffprim  + dmeffsec  )
flopy.mf6.ModflowGwtdsp(
    gwt,
    alh      = 0,
    ath1     = 0,
    diffc    = equivdiffc,
    xt3d_off = False,
    filename = "{}.dsp".format(gwtname),
)

# mst
flopy.mf6.ModflowGwtmst(
    gwt,
    porosity          = porosity,
    first_order_decay = False,
    decay             = None,
    decay_sorbed      = None,
    sorption          = None,
    bulk_density      = None,
    distcoef          = None,
    filename="{}.mst".format(gwtname),
)

# ssm
sourcerecarray = None 
flopy.mf6.ModflowGwtssm(
    gwt,
    sources=sourcerecarray,
    filename="{}.ssm".format(gwtname),
)

# gwt
concentration_file = "{}.ucn".format(gwtname)
if SAVEALL:
    saverecord=[ ('CONCENTRATION', 'ALL') ]
else:
    saverecord=[ ('CONCENTRATION', 'FIRST'), ('CONCENTRATION','LAST') ]
flopy.mf6.ModflowGwtoc(
    gwt,
    concentration_filerecord=concentration_file,
    concentrationprintrecord=[
        ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
    ],
    saverecord=saverecord, 
    filename="{}.oc".format(gwtname),
)

# gwfgwt
flopy.mf6.ModflowGwfgwt(
    sim,
    exgtype="GWF6-GWT6",
    exgmnamea=base_model_name,
    exgmnameb=gwtname,
    filename="{}.gwfgwt".format(base_model_name),
)

# write sim 
write_simulation = args.write
if write_simulation:
    sim.write_simulation()

# init modflow-api 
libdir = config.FOLDERS['lib']
mf6    = xmipy.XmiWrapper(
        os.path.join( libdir, 'libmf6.so' ),  # linux 
        #os.path.join( libdir, 'libmf6.dll' ), # wos
        working_directory = sim_ws
    )

# run sim 
if args.run:
    gwtname = gwtname.upper()
    mf6.set_int("ISTDOUTTOFILE", 0)
    mf6.initialize( os.path.join( sim_ws, 'mfsim.nam' ) )
    current_time = mf6.get_current_time()
    end_time     = mf6.get_end_time()
    timecounter  = 0
    while current_time < end_time:

        # UPDATE
        mf6.update()

        # GET CONCENTRATION
        scalar = mf6.get_value(gwtname+'/X')

        # Calculate nonlinear function
        nlfunc = scalar/np.sqrt(scalar**2+1)

        # Calculate equivalent diffusion
        equivdiffc  = 0.5*( ( dmeffprim  - dmeffsec  )*nlfunc + dmeffprim  + dmeffsec  )

        # Update equivalent diffusion
        mf6.set_value( gwtname+'/DSP/DIFFC' , equivdiffc  )

        # NEXT TIME
        current_time = mf6.get_current_time()

    # FINALIZE
    mf6.finalize()

