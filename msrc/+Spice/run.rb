#!/usr/bin/env ruby

require 'progressbar'

netlist = 'nmos.cir'
netlist_tmp = "#{ netlist }_tmp"

if ARGV.count > 2
  netlist_out = ARGV[2]
else
  netlist_out = File.basename(netlist, '.cir') + '.out'
end

dir = File.dirname(__FILE__)

netlist = File.open(File.join dir, netlist).read

#
# Channel length
#

if ARGV.count > 0
  Lcount = ARGV[0].to_i
else
  Lcount = 100
end

Lnom = 22e-9
Ldev = 0.05 * Lnom

Lmin = Lnom - Ldev
Lrange = 2 * Ldev

if Lcount == 1
  L = [ Lnom ]
else
  L = Array.new(Lcount)

  (0...Lcount).to_a.each do |i|
    L[i] = "%.4e" % (Lmin + i * Lrange / (Lcount - 1))
  end
end

#
# Temperature
#

if ARGV.count > 1
  Tcount = ARGV[1].to_i
else
  Tcount = 100
end

Tnom = 27

Tmin = Tnom
Trange = 100 - Tmin

if Tcount == 1
  T = [ Tnom ]
else
  T = Array.new(Tcount)

  (0...Tcount).to_a.each do |i|
    T[i] = "%.2f" % (Tmin + i * Trange / (Tcount - 1))
  end
end

#
# Simulation
#

output = File.open(netlist_out, 'w')

output.puts "L\tT\tI"

pb = ProgressBar.new('Ngspice', L.length * T.length);

L.each do |l|
  T.each do |t|
    netlist0 = netlist.gsub(/L\s*=\s*[^\s]+/, "L = #{ l }")
    netlist0 = netlist0.gsub(/Temp\s*=\s*[^\s]+/, "Temp = #{ t }")

    File.open(File.join(dir, netlist_tmp), 'w') { |f| f.write netlist0 }
    pipe = IO.popen("cd #{ dir } && ngspice -b #{ netlist_tmp } 2> /dev/null")
    lines = pipe.readlines
    pipe.close

    lines.each do |line|
      next unless line.match /@m1\[id\]\s*=\s*(.*)/
      output.puts "#{ l }\t#{ t }\t#{ $1 }"
      break
    end

    pb.inc
  end
end

puts

output.close

File.delete File.join(dir, netlist_tmp) rescue true
File.delete File.join(dir, 'bsim4.out') rescue true
