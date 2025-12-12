export Parameter, optimize

##########################################################################################
#   Parameter
##########################################################################################
abstract type Parameter end

##########################################################################################
#   Implementation: optimize
##########################################################################################
function optimize(param::Parameter, qm::QuantificationMeasure, args...)
    error("The 'optimize' is not implemented to the parameter '$(typeof(param))' for the measure '$(typeof(qm))'.")
end