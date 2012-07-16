#!/usr/bin/env ruby

require 'progressbar'

netlist = 'nmos.cir'
netlist_tmp = "#{ netlist }_tmp"
netlist_out = File.basename(netlist, '.cir') + '.out'

dir = File.dirname(__FILE__)

netlist = File.open(File.join dir, netlist).read

params = {
  'Temp' => (27..100)
}

output = File.open(netlist_out, 'w')

params.keys.each do |param|
  pb = ProgressBar.new(param, params[param].count)

  output.puts "#{ param }\tValue"

  params[param].each do |value|
    netlist0 = netlist.gsub(/#{ param }\s*=\s*[^\s]+/, "#{ param } = #{ value }")
    File.open(File.join(dir, netlist_tmp), 'w') { |f| f.write netlist0 }
    pipe = IO.popen("cd #{ dir } && ngspice -b #{ netlist_tmp } 2> /dev/null")
    lines = pipe.readlines
    lines.each do |line|
      next unless line.match /@m1\[id\]\s*=\s*(.*)/
      output.puts "#{ value }\t#{ $1 }"
      break
    end
    pb.inc
  end

  output.puts
  puts
end

output.close

File.delete File.join(dir, netlist_tmp)
File.delete File.join(dir, 'bsim4.out')
