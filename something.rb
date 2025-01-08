# name this Ruby file something.rb
# execute by running this in a terminal:
#   ruby something.rb
a,b = 1,1
puts a
puts b
10.times do
  puts a+b
  a,b = b,a+b
end