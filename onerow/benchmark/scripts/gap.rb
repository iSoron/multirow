#!/usr/bin/env ruby
require 'yaml'
yaml = YAML::load(File.open(ARGV[0]))

int_format = "%d,"
number_format = "%.2f,"
string_format = "%s,"

presolve_opt = yaml['sol_value'][0] 
mir_opt = yaml['sol_value'][1]
wedges_opt = yaml['sol_value'][2]
mip_opt = yaml['mip_value']

if presolve_opt == mip_opt then
    print (string_format * 6 + "\n") % (["---"] * 6)
else

    orig_gap = (mip_opt-presolve_opt) / mip_opt.abs
    mip_perf  = (mir_opt-presolve_opt) / (mip_opt-presolve_opt)
    w_perf    = (wedges_opt-presolve_opt) / (mip_opt-presolve_opt)

    n_cuts_mir = yaml['n_added_cuts']['depth'][0]
    n_cuts_w   = yaml['n_added_cuts']['depth'].values.reduce(:+) - n_cuts_mir

    print int_format % n_cuts_mir
    print int_format % n_cuts_w
    print number_format % (100 * orig_gap)
    print number_format % (100 * mip_perf)
    print number_format % (100 * w_perf)
    
    if  presolve_opt == wedges_opt
        print (string_format * 2) % (["---"] * 2)
    else
        mir_contrib = (mir_opt-presolve_opt) / (wedges_opt-presolve_opt)
        w_contrib   = (wedges_opt-mir_opt) / (wedges_opt-presolve_opt)

        print number_format % (100 * mir_contrib)
        print number_format % (100 * w_contrib)
    end

    w_improv = (wedges_opt-mir_opt) / (mip_opt-mir_opt)
    print (number_format ) % (100 * w_improv.abs)
    print (int_format + "\n") % yaml['timers'][2]
end
