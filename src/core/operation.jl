export Operation, operate

##########################################################################################
#   Operations with microstates
##########################################################################################
abstract type Operation end

##########################################################################################
#   Function to operate
##########################################################################################
function operate(op::Operation)
    error("Not implemented operation without arguments: $(typeof(op))")
end