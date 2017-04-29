#!/usr/bin/env ruby
require 'yaml'
yaml = YAML::load(File.open(ARGV[0]))

int_format = "%d,"
number_format = "%.2f,"
string_format = "%s,"
blank = "---,"

print int_format % yaml['n_generated_cuts']['round'][1] rescue print blank
print int_format % yaml['n_generated_cuts']['round'][2] rescue print blank
print number_format % (yaml['cut_speed']['round_1'] * 1000) rescue print blank
print number_format % (yaml['cut_speed']['round_2'] * 1000) rescue print blank
print number_format % yaml['trivial_lifting']['average_m'] rescue print blank
print "\n"
