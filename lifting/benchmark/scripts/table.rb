#!/usr/bin/env ruby
# encoding: UTF-8
require 'yaml'
require 'bigdecimal'

def sum(a)
    a.inject(0){ |accum, i| accum + i }
end

def average(a)
    sum(a) / a.length.to_f
end

def gmean(a)
    prod = BigDecimal.new 2
    a.each { |v| prod *= (BigDecimal.new(v,5))}
    prod ** (1.0 / a.size)
end

def median(a)
    a.sort[a.length / 2]
end

def standard_deviation(a)
    m = average(a)
    sum = a.inject(0){ |accum, i| accum + (i - m) ** 2 }
    Math.sqrt(sum / (a.length - 1).to_f)
end

files = []
filenames = []
ARGV.each_with_index do |filename,idx|
    filenames[idx] = filename
            .gsub(/.*\//,"")
            .gsub(".yaml", "")
            .gsub(/[^A-Za-z0-9]/, "")
    files[idx] = YAML::load(File.open(filename))
end

fail_avg = []
fail_count = []
success_count = []
times = []
all_times = []
ratios_to_best = []

n_instances = files[0]['cpu_time'].length

files.each_with_index do |f,idx|

    fail_count[idx] = 0
    success_count[idx] = 0
    times[idx] = []
    all_times[idx] = []
    ratios_to_best[idx] = []

    for i in 0..(n_instances-1) do
        time = f['cpu_time'][i]
        next if time.nil?
        times[idx].push(time)
        all_times[idx].push(time)

        wrong = f['wrong_answer'][i]
        if(wrong == 1)
            fail_count[idx] += 1
        else
            success_count[idx] += 1
        end
    end

    fail_avg[idx] = fail_count[idx].to_f / n_instances
end

best_percentage = [0] * files.length

BIG_M = 1000000000
for i in 0..(n_instances-1) do
    best_time = BIG_M
    files.each_with_index do |f,idx|
        time = all_times[idx][i] || BIG_M
        best_time = [best_time, time].min
    end

    files.each_with_index do |f,idx|
        time = times[idx][i] || BIG_M
        best_percentage[idx] += 1 if time == best_time

        next if best_time <= 0
        ratios_to_best[idx].push(time / best_time)
    end
end

best_percentage = best_percentage.map{|b| b.to_f / n_instances}

def print_line(title, format, arr)
    print "%s," % title
    arr.each do |a|
        print format % a
    end
    print "\n"
end

print_line("metric", "%s,", filenames)
print_line("Average (ms)", "%.3f,", times.map{|t| average(t)})
print_line("Median (ms)", "%.3f,", times.map{|t| median(t)})
print_line("Maximum (ms)", "%.3f,", times.map{|t| t.max})
print_line("Failure Rate", "%.1f \\%%,", fail_avg.map{|x| 100.0 * x})
print_line("Best", "%.1f \\%%,", best_percentage.map{|x| 100 * x})
print_line("Avg Ratio to Best", "%.3f,", ratios_to_best.map{|t| average(t)})
