#!/usr/bin/env python

#################################################################################
#  RELIES UPON CDAT MODULES 
#  PRODUCE SIMPLE PORTRAIT PLOT BASED ON RESULTS FROM PCMDI'S METRICS PACKAGE 
#  WRITTEN BY P. GLECKLER, PCMDI/LLNL
#  FURTHER TWEAKED BY C. DOUTRIAUX
#  ISSUE TO BE RESOLVED WITH UV-CDAT PLANNED IMPROVEMENTS: ABILITY TO SAVE PLOTS IN MULTIPLE FORMATS  
#################################################################################

# CHECK TO VERIFY VCS MODULE IS INSTALLED
try:
    import vcs
except:
    raise RuntimeError("Package: vcs, not found to generate portrait plots")

# STANDARD PYTHON MODULES
import glob
import json
import os
import sys
import re
try:
    import numpy as np
except:
    raise RuntimeError("Package: numpy, not found.")

# CDAT MODULES
try:
    import pcmdi_metrics.graphics.portraits
except:
    raise RuntimeError("Package: pcmdi_metrics.graphics.portraits not found.")    
try:
    import MV2
except:
    raise RuntimeError("Package: MV2 not found.")
import argparse
from pprint import pprint

def create_plot(args):
    # CREATE VCS OBJECT AS A PORTAIT PLOT AND LOAD PLOT SETTINGS FOR TEST CASE 
    x=vcs.init()
    x.portrait()

    #######################
    # PORTRAIT PLOT PARAMETER SETTINGS
    P=pcmdi_metrics.graphics.portraits.Portrait()
    # Turn off verbosity
    P.verbose = False

    P.PLOT_SETTINGS.levels = [-1.e20,-.5,-.4,-.3,-.2,-.1,0.,.1,.2,.3,.4,.5,1.e20]
    P.PLOT_SETTINGS.x1 = .1
    P.PLOT_SETTINGS.x2 = .85
    P.PLOT_SETTINGS.y1 = .22
    P.PLOT_SETTINGS.y2 = .95
    #P.PLOT_SETTINGS.xtic2.y1=P.PLOT_SETTINGS.y1
    #P.PLOT_SETTINGS.xtic2.y2=P.PLOT_SETTINGS.y2
    #P.PLOT_SETTINGS.ytic2.x1=P.PLOT_SETTINGS.x1
    #P.PLOT_SETTINGS.ytic2.x2=P.PLOT_SETTINGS.x2
    # Logo is simply a vcs text object, set to None for off
    P.PLOT_SETTINGS.logo = None
    # timestamp is simply a vcs text object, set to None for off
    P.PLOT_SETTINGS.time_stamp = None
    P.PLOT_SETTINGS.draw_mesh='n'
    x.scriptrun(os.path.join(sys.prefix,"share","pmp","graphics","vcs","portraits.scr"))
    P.PLOT_SETTINGS.colormap = 'bl_rd_12'
    cols=vcs.getcolors(P.PLOT_SETTINGS.levels,range(144,156),split=1)
    P.PLOT_SETTINGS.fillareacolors = cols
    P.PLOT_SETTINGS.parametertable.expansion = 100 

    # LIST OF VARIABLES TO BE USED IN PORTRAIT PLOT
    # REDUCED LIST FOR TEST CASE
#    vars = ['pr','tas','psl','rlut','rsut','ua-850','ua-200','va-850','va-200','zg-500']
    variables = args.variables.split(",")

    # LOAD METRICS DICTIONARIES FROM JSON FILES FOR EACH VAR AND STORE AS A SINGLE DICTIONARY
    var_cmip5_dics = {}
    mods = set()

    # CMIP5 METRICS RESULTS - CURRENTLY USING FOR CONTROL SIMULATIONS
    if args.cmip_compare:
        cmip_json_files = glob.glob(os.path.join(args.cmipmetrics,'CMIP_metrics_results','CMIP5',args.mode,'*json'))

    # ADD GFDL JSON FILES... 
    # This is pretty hard coded might want to consider more magic
    json_files = []
    input_experiments = []
    for input_dir in args.argv:
        dir_files = glob.glob(os.path.join(input_dir,'*.json'))
        with open(dir_files[0], "r") as f:
            data = json.load(f)

        input_experiments.append(data["RESULTS"].keys()[0])
        json_files += dir_files
        

    if not json_files:
        sys.stderr.write("ERROR: No JSON files found in the input directory.\n")
        sys.exit(1)

    json_files += cmip_json_files
    #Aparna edits, to work in jupyter
    EXCLUDE_MODELS = list()
    if args.exclude:
        EXCLUDE_MODELS = [ model.lower() for model in args.exclude.split(",") ]
    # CONSTRUCT PYTHON DICTIONARY WITH RESULTS METRICS USED IN PORTRAIT  
    non_mods = ["GridInfo","References","RegionalMasking","metrics_git_sha1","uvcdat_version"]
    for fnm in json_files:
        f=open(fnm)
        d = json.load(f)
        var = os.path.basename(fnm).split("_")[0]
        for m in d['RESULTS'].keys():
            # Skip non model bits`
            if m not in non_mods:
                if args.gfdl_compare:
                    if 'GFDL' in m or m in input_experiments:
                        if m.lower() not in EXCLUDE_MODELS:
                            mods.add(m)
                        else:
                            del(d["RESULTS"][m])
                    else:
                        del(d["RESULTS"][m])
                elif m.lower() in EXCLUDE_MODELS:
                    del(d["RESULTS"][m])
                else:
                    mods.add(m)
                    
            else:
                # Let's clean it up
                del(d[m])
            if var_cmip5_dics.has_key(var):
                var_cmip5_dics[var]['RESULTS'].update(d["RESULTS"])
            else:
                var_cmip5_dics[var]=d
    
    mods = sorted([ mod for mod in mods if mod not in input_experiments ])
#    gfdl_model = mods.pop(0)
    mods = input_experiments + mods

#    mods.insert(0,gfdl_model)

    # ORGANIZE METRICS INTO A VARIABLES X MODELS MATRIX 
    out1_rel,out2_rel,out3_rel,out4_rel = [np.ma.masked_all((len(variables),len(mods)),np.float32) for _ in range(4)] ; # Define arrays to fill

    # LOOP OVER VARIABLE
    for vn, var in enumerate(variables):
        vals1,vals2,vals3,vals4 = [[] for _ in range(4)]
        # LOOP OVER MODEL
        for mn,mod in enumerate(mods):
            try:
                out1_rel[vn,mn] = float(var_cmip5_dics[var]['RESULTS'][mod]["defaultReference"]["r1i1p1"]["global"]['rms_xy_djf'])
                out2_rel[vn,mn] = float(var_cmip5_dics[var]['RESULTS'][mod]["defaultReference"]["r1i1p1"]["global"]['rms_xy_mam'])
                out3_rel[vn,mn] = float(var_cmip5_dics[var]['RESULTS'][mod]["defaultReference"]["r1i1p1"]["global"]['rms_xy_jja'])
                out4_rel[vn,mn] = float(var_cmip5_dics[var]['RESULTS'][mod]["defaultReference"]["r1i1p1"]["global"]['rms_xy_son'])
            except Exception,err:
                pass
    
        med_rms1 = np.ma.median(out1_rel[vn,:])
        med_rms2 = np.ma.median(out2_rel[vn,:])
        med_rms3 = np.ma.median(out3_rel[vn,:])
        med_rms4 = np.ma.median(out4_rel[vn,:])
   
        out1_rel[vn,:]=(out1_rel[vn,:]-med_rms1)/med_rms1
        out2_rel[vn,:]=(out2_rel[vn,:]-med_rms2)/med_rms2
        out3_rel[vn,:]=(out3_rel[vn,:]-med_rms3)/med_rms3
        out4_rel[vn,:]=(out4_rel[vn,:]-med_rms4)/med_rms4

    # ADD SPACES FOR LABELS TO ALIGN AXIS LABELS WITH PLOT
    yax = [ m.encode('utf-8')+" " for m in mods ]
    xax = [ v+" " for v in variables ]
  
    # Convert to MV so we can decorate
    out1_rel = MV2.array(out1_rel)
    out2_rel = MV2.array(out2_rel)
    out3_rel = MV2.array(out3_rel)
    out4_rel = MV2.array(out4_rel)

    # GENERATE PLOT 
    P.decorate(out1_rel,xax,yax)
    P.decorate(out2_rel,xax,yax)
    P.decorate(out3_rel,xax,yax)
    P.decorate(out4_rel,xax,yax)

    # PLOT
    P.plot(out1_rel,x=x,multiple=1.4)
    P.plot(out2_rel,x=x,multiple=2.4)
    P.plot(out3_rel,x=x,multiple=3.4)
    P.plot(out4_rel,x=x,multiple=4.4)
    # END OF PLOTTING

    # SAVE PLOT
    fileName = ".".join([args.descriptor, args.format.lower()])
    if re.match("png", args.format, re.IGNORECASE):
        x.png(os.path.join(args.output,fileName))
    elif re.match("ps", args.format, re.IGNORECASE):
        x.postscript(os.path.join(args.output,fileName), textAsPaths=args.compress)
    elif re.match("pdf", args.format, re.IGNORECASE):
        x.pdf(os.path.join(args.output,fileName),textAsPaths=args.compress)
    elif re.match("eps", args.format, re.IGNORECASE):
        x.eps(os.path.join(args.output,fileName),textAsPaths=args.compress)
    else:
        sys.stderr.write("ERROR: I don't know that format: %s\n" % args.format)
        sys.exit(1)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Portrait plot generation tool.")
#    parser.add_argument("-i", "--input", help="The input directory to the PCMDI Metrics.", required=True)
    parser.add_argument("-d", "--descriptor", help="The experiment name.", required=True)
    parser.add_argument("-o", "--output", help="The output directory for the portrai plot.",
                        default=os.getcwd(), required=False)
    parser.add_argument("-f", "--format", help="The output file type for the protrait plot (ps, pdf, png)",
                        default="png", required=False)
    parser.add_argument("-C", "--compress", action="store_false", help="Compress the output.",
                        required=False)
    parser.add_argument("-g", "--gfdl_compare", help="Option to only compare against GFDL CMIP models.", action="store_true", default=False)
    parser.add_argument("-c", "--cmip_compare", help="Option to only compare against the cmip_models.", action="store_false", default=True) 
    parser.add_argument("-v", "--variables", help="Comma delimited list of variables to show.",
                        default="pr,psl,tas,rlut,rsut,ua-850,ua-200,va-850,va-200,zg-500", required=False)
    parser.add_argument("-x", "--exclude", help="Comma delimited list to exclude CMIP5 models.", required=False)
    parser.add_argument("-m", "--mode", help="The option to change which version of the data to compare against.", default="amip")
<<<<<<< HEAD
=======
    parser.add_argument("-s", "--cmipMetrics", help="specifies metrics location  computed from other CMIP models. specify path till share directory only", default="/home/a1r/metrics_latest_v1p1p1/share/")
>>>>>>> 76717f2ec0e6c273e9b68e60eaf5a8f0057abe64
    parser.add_argument("argv", nargs=argparse.REMAINDER)
    args = parser.parse_args()
 
    if len(args.argv) == 0:
        raise RuntimeError("Sorry, but it appears that you've not specified any input directories!")
    
    EXCLUDE_MODELS = list()
    if args.exclude:
        EXCLUDE_MODELS = [ model.lower() for model in args.exclude.split(",") ]

    GFDL_COMPARE = args.gfdl_compare
    CMIP_COMPARE = args.cmip_compare
    if len(args.argv) == 1:
        # only one model provided on the commandline, so lets compare against CMIP!
        CMIP_COMPARE = True

    if GFDL_COMPARE and not CMIP_COMPARE:
        CMIP_COMPARE = True

    for input_dir in args.argv:
        if not os.path.exists(input_dir):
            sys.stderr.write("ERROR: The input directory for the PCMDI Metrics JSON files doesn't exist.\n")
            sys.exit(1)

    try:
        if not os.path.exists(args.output):
            os.mkdir(args.output)
    except:
        sys.stderr.write("ERROR: Unable to make the output directory for the portrait plots.\n")
        sys.exit(1)

    create_plot(args)
