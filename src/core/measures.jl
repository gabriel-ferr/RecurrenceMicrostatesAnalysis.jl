export QuantificationMeasure, measure

##########################################################################################
#   Type: Quantification Measure
##########################################################################################
abstract type QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################

function measure(ms::QuantificationMeasure, agrs...)
    error("There isn't a 'measure' implementation for '$(typeof(ms))'.")
end