export RMACore, histogram, distribution

##########################################################################################
#   RMACore
##########################################################################################
abstract type RMACore end

##########################################################################################
#   Implementations: histogram & distribution
##########################################################################################
function histogram(
    core::C,
    x,
    y;
    threads = Threads.nthreads()
) where { C <: RMACore }
    throw("The RMA core of type '$(typeof(core))' is not implemented for input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end

function distribution(
    core::C,
    x,
    y;
    threads = Threads.nthreads()
) where { C <: RMACore }
    throw("The RMA core of type '$(typeof(core))' is not implemented for input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end

