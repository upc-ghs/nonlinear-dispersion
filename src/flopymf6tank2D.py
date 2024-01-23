'''
FloPy routine writing config files for the simulation of reactive 
transport of a binary chemical system in a two-dimensional 
porous media tank. 

Scenario considers a system initially saturated with the 
primary species (c_p), and the secondary (c_s) is injected at the 
center of the system.

Problem is solved in terms of a component u = c_p - c_s, 
with an equivalent nonlinear dispersion. 
'''

# Dependencies
import os
import sys
import numpy as np
import flopy
import matplotlib.pyplot as plt
import argparse
import xmipy
import scipy.special as spc

# self
import config
import dispersion

#############
# Arguments #
#############
parser = argparse.ArgumentParser( description='Setup simulations for experiment.' )
parser.add_argument( '--write'    , action='store_true', help='write simulation files.' )
parser.add_argument( '--run'      , action='store_true', help='executes simulation, requires mf6lib.' )
parser.add_argument( '--force'    , action='store_true', help='skip the prompt after verification of time step.' )
parser.add_argument( '--cprofile' , action='store_true', help='displays a transverse profile of concentration.' )
args   = parser.parse_args()


#########################
# General configuration #
#########################
HIGHRES        = False # adopt high resolution parameters
targetvelocity = 10 # m/day 

# total time in hours 
if targetvelocity == 0.5: 
    tlength       = 20.0 
if targetvelocity == 1: 
    tlength       = 10.0   
if targetvelocity == 10: 
    tlength       = 1.0   
if targetvelocity == 20: 
    tlength       = 0.5
if HIGHRES: 
    tn_time_steps = 15000
else:
    tn_time_steps = 1500


stressperiods = [
        {
            'id': 0,
            'length'       : tlength,
            'n_time_steps' : tn_time_steps, 
            'ts_multiplier': 1,
            'steady_state' : True,
            'printfreqts'  : 1, 
        },
    ]
experiment_folder = config.FOLDERS['sims']
base_model_name   = 'mf6model'

# component
component         = 0.5
dpdsdiff          = 0.1 # dp/ds
sim_name          = 'NLCOMP1D.dpds.'+str(dpdsdiff)+'.u'+str(component)
davg              = 1e-9*1e4*3600   # avg diffusion, the same for all sims
dmscale           = 1e-8*1e4*3600   # scale for time step selection 
daqsec            = 2*davg/( 1 + dpdsdiff )
daqprim           = dpdsdiff*daqsec
porosity          = 0.35
delta             = 6
betaT             = 0.5
ck                = 0.01             # np.sqrt(ms) ( Hazen, 1892 )
diam              = 0.0015           # meters (1.5mm)
hk                = diam**2/(ck**2)  # ( Hazen, 1892 )
hk                = hk*100*3600      # m/s
sim_ws            = os.path.join( experiment_folder, sim_name )

print( ' - INFO:  Grain diameter  : ', diam*1000  , 'mm   ' ) 
print( ' - INFO:  Hydraulic cond  : ', hk         , 'm/s  ' )


# Init sim and configure stress periods
sim = flopy.mf6.MFSimulation(
        sim_name=sim_name,
        sim_ws=sim_ws,
        exe_name='mf6'
    )
# tdis package
perioddata = []
for sp in stressperiods:
    perioddata.append( [sp['length'], sp['n_time_steps'], sp['ts_multiplier']] )
tdis = flopy.mf6.ModflowTdis(
    sim,
    nper=len(stressperiods),
    perioddata=perioddata,
    time_units='hours',
)

# IMS
hclose = 1e-5
ninner = 100
nouter = 50
relax  = 0.97
complexity = 'COMPLEX'
ims = flopy.mf6.ModflowIms(
        sim,
        pname               = 'ims',
        print_option        = 'SUMMARY',
        complexity          = complexity,
        outer_dvclose       = hclose, 
        outer_maximum       = nouter, 
        inner_maximum       = ninner, 
        inner_dvclose       = hclose, 
        linear_acceleration = 'BICGSTAB',
    )

# gwf
gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=base_model_name,
        save_flows=True
    )
# dis package 
LTANK  = 40.0  # cm
nlay   = 1
if HIGHRES:
    wcells = 31    # number of cells injecting
    nrow   = 250
    ncol   = 1250
else:
    wcells = 12    # number of cells injecting
    nrow   = 160
    ncol   = 800
delc   = LTANK/ncol
delz   = delc
delr   = delc
dis    = flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay, nrow=nrow, ncol=ncol,
            delr=delr, delc=delc, top=delz, botm=0,
            length_units='CENTIMETERS',
        )

# npf package
npf = flopy.mf6.ModflowGwfnpf(
    gwf,
    save_specific_discharge=True,
    icelltype=0,
    k=hk,
)

# sto package
sto = flopy.mf6.ModflowGwfsto(gwf, ss=0,sy=0, steady_state=True)

# chd 
# Darcy's law
print( ' - INFO:  Target velocity : ', targetvelocity, '  m/d  ')
targetvelocity  = targetvelocity*100/24   # cm/hr
print( ' - INFO:  Target velocity : ', targetvelocity, ' cm/hr ')
qflow           = targetvelocity*porosity # cm/hr
gradh           = qflow/hk                # deltah/length
deltah          = gradh*LTANK
href            = 10 # could be anything 
hdin            = href + deltah
hdout           = href
totalflowrate   = qflow*delc**2*nrow
print( ' - INFO:  Head gradient   : ', gradh  , ' -  ')
print( ' - INFO:  Delta head      : ', deltah , ' cm ')
print( ' - INFO:  Head inlet      : ', hdin   , ' cm ')
print( ' - INFO:  Head outlet     : ', hdout  , ' cm ')
print( ' - INFO:  Flow rate       : ', totalflowrate , ' cm3/hr ' )
print( ' - INFO:  Pore volume     : ', porosity*nrow*ncol*nlay*delc**3/totalflowrate    , ' hr  ')
print( ' - INFO:  Pore volume     : ', porosity*nrow*ncol*nlay*delc**3/totalflowrate*60 , ' min ')

# assign chd for in flow cells 
# center cells inject secondary species
# the rest inject the primary species
chdspd     = []
thesecells = [ int(0.5*(nrow - wcells)) + iw for iw in range(wcells) ]
for ir in range(nrow):
    if ir in thesecells:
        chdspd.append( [ (0, ir, 0), hdin, -1*component  ] ) # inject secondary at the center
    else:
        chdspd.append( [ (0, ir, 0), hdin, component ] )     # inject primary in the rest
chdin = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        maxbound=len(chdspd),
        auxiliary="CONCENTRATION",
        pname='CHD-1'
    )
# assign chd for outflow cells 
chdspd     = []
for ir in range(nrow):
    chdspd.append( [ (0, ir, ncol-1), hdout ] )
chdout = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        maxbound=len(chdspd),
        pname='CHD-2'
    )

# ic (linear gradient)
icst = np.zeros(shape=(nlay,nrow,ncol))
xarr = np.arange(0.5*delc,LTANK+0.5*delc,delc)
for ir in range(nrow):
    icst[0,ir,:] = -1*gradh*xarr + hdin
flopy.mf6.ModflowGwfic(
        gwf,
        strt=icst,
    )

# oc
head_file   = base_model_name + '.hds'
budget_file = base_model_name + '.bud'   
nsteps      = stressperiods[0]['n_time_steps']
tlength     = stressperiods[0]['length']
deltat      = tlength/nsteps
printfreqts = stressperiods[0]['printfreqts']
times       = np.arange(0,tlength+deltat,deltat)
printtimes  = times[::printfreqts]
saverecord  = [ ('HEAD','LAST'), ('BUDGET','LAST')]
flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord  =head_file,
    budget_filerecord=budget_file,
    saverecord=saverecord
)

#######
# GWT #
#######
gwtname = "gwt-" + base_model_name
gwt     = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file="{}.nam".format(gwtname),
    )
gwt.name_file.save_flows = True

# create iterative model solution and register the gwt model with it
imsgwt = flopy.mf6.ModflowIms(
    sim,
    print_option        = "SUMMARY",
    outer_dvclose       = hclose,
    outer_maximum       = nouter,
    under_relaxation    = "NONE",
    inner_maximum       = ninner,
    inner_dvclose       = hclose,
    linear_acceleration = "BICGSTAB",
    complexity          = complexity,
    filename            = "{}.ims".format(gwtname),
)
sim.register_ims_package(imsgwt, [gwt.name])

# disgwt
disgwt = flopy.mf6.ModflowGwfdis(
        gwt,
        nlay=nlay, nrow=nrow, ncol=ncol,
        delr=delr, delc=delc, top=delz, botm=0,
        length_units='CENTIMETERS',
    )

# adv
scheme = "TVD"
flopy.mf6.ModflowGwtadv(
        gwt,
        scheme=scheme,
        filename="{}.adv".format(gwtname)
    )

#ic conc
#saturated with primary species
sconc = component
flopy.mf6.ModflowGwtic(
    gwt,
    strt=sconc,
    filename="{}.ic".format(gwtname)
)

# some parameters for dispersion
deltax      = delc 
difffactor  = daqprim/daqsec
tscaledaqp  = deltax**2/(dmscale)    # cm2/hr
tscaledaqs  = deltax**2/(dmscale)    # cm2/hr
tscaleadv   = deltax/targetvelocity  # hr 
ct          = 0.1
displ       = 0.5*diam*100*targetvelocity # cm2/hr
displphip   = porosity*daqprim + 0.5*diam*100*targetvelocity # cm2/hr
displphis   = porosity*daqsec  + 0.5*diam*100*targetvelocity # cm2/hr
tscaledlp   = deltax**2/(displphip)    # cm2/hr
tscaledls   = deltax**2/(displphis)    # cm2/hr
trec        = np.min([ct*np.mean(tscaledlp), ct*np.mean(tscaledls), ct*tscaleadv])
pecletprim  = targetvelocity*diam*100/daqprim
pecletsec   = targetvelocity*diam*100/daqsec

# equivalent dispersion parameters
# primary species
alphalprim  = 0.5*diam*100 # cm
alphatprim  = (diam*100/pecletprim)*dispersion.nonlinear_dispersion_function( pecletprim, delta, betaT )
dmeffprim   = porosity*daqprim
# secondary species
alphalsec   = 0.5*diam*100 # cm
alphatsec   = (diam*100/pecletsec )*dispersion.nonlinear_dispersion_function( pecletsec , delta, betaT )
dmeffsec    = porosity*daqsec

# dsp package
# initialize equivalent dispersion for the component
nlfunc      = sconc/np.sqrt(sconc**2+1) # u/sqrt(u**2+1)
equivdiffc  = 0.5*( ( dmeffprim  - dmeffsec  )*nlfunc + dmeffprim  + dmeffsec  )
equivalphal = 0.5*( ( alphalprim - alphalsec )*nlfunc + alphalprim + alphalsec )
equivalphat = 0.5*( ( alphatprim - alphatsec )*nlfunc + alphatprim + alphatsec )
flopy.mf6.ModflowGwtdsp(
    gwt,
    alh      = equivalphal,
    ath1     = equivalphat,
    diffc    = equivdiffc,
    xt3d_off = False,
    filename = "{}.dsp".format(gwtname),
)

print( ' # Dispersion analysis:    ' )
print( '    - deltat (given)      : ', stressperiods[0]['length']/stressperiods[0]['n_time_steps'] )
print( '    - deltat (recommended): ', trec )
print( '    - nsteps (given)      : ', stressperiods[0]['n_time_steps'] ) 
print( '    - nsteps (recommended): ', stressperiods[0]['length']/trec ) 
print( '    - effective CT        : ', (stressperiods[0]['length']/stressperiods[0]['n_time_steps'])/trec*ct  ) 
print( 'Continue ... ? (press c to continue, or q to exit)' )
if not args.force:
    import pdb
    pdb.set_trace()


# mst package
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
sourcerecarray = [('CHD-1','AUX','CONCENTRATION')]
flopy.mf6.ModflowGwtssm(
    gwt,
    sources=sourcerecarray,
    filename="{}.ssm".format(gwtname),
)

# gwt oc 
concentration_file = "{}.ucn".format(gwtname)
saverecord         = [ ('CONCENTRATION', 'FIRST'), ('CONCENTRATION','LAST') ]
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

#############
# Write MF6 #
#############
write_simulation = args.write
if write_simulation:
    sim.write_simulation()

# MODFLOW API
libdir = config.FOLDERS['lib']
mf6    = xmipy.XmiWrapper(
        os.path.join( libdir, 'libmf6.so' ), 
        working_directory = sim_ws
    )

###############
# Execute MF6 #
###############
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

        # Calculate equivalent diffusion and dispersivities 
        equivdiffc  = 0.5*( ( dmeffprim  - dmeffsec  )*nlfunc + dmeffprim  + dmeffsec  )
        equivalphal = 0.5*( ( alphalprim - alphalsec )*nlfunc + alphalprim + alphalsec )
        equivalphat = 0.5*( ( alphatprim - alphatsec )*nlfunc + alphatprim + alphatsec )

        # Update equivalent dispersivities and diffusion
        mf6.set_value( gwtname+'/DSP/DIFFC' , equivdiffc  )
        mf6.set_value( gwtname+'/DSP/ALH'   , equivalphal )
        mf6.set_value( gwtname+'/DSP/ATH1'  , equivalphat )
        mf6.set_value( gwtname+'/DSP/ATH2'  , equivalphat )

        # NEXT TIME
        current_time = mf6.get_current_time()

    # FINALIZE
    mf6.finalize()


if args.cprofile:

    import matplotlib.colors as colors
    import matplotlib.cm as cm 
    
    gwtname  = gwtname.lower()
    layer_id = 0
    fname    = os.path.join( os.getcwd(), sim_ws, gwtname+'.ucn')
    ucnobj   = flopy.utils.HeadFile(fname, text='CONCENTRATION', precision='double')
    data     = ucnobj.get_alldata()
    xindex   = -50 

    # Figure handling
    modelgrid      = gwf.modelgrid
    domainwidth    = nrow*delr
    deltay         = delr
    yvector        = modelgrid.xycenters[1] - domainwidth*0.5 + deltay
    plt.plot( yvector, data[-1,0,:,xindex], color='r' )
    ax = plt.gca()
    ax.set_xlabel('y')
    ax.set_ylabel('C')

    plt.savefig(os.path.join(sim_ws,'cprof'+sim_name+'.png'))
