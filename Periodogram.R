# Author: Patrick Reidy
# Email:  <patrick.francis.reidy@gmail.com>





################################################################################
# Class definition                                           Periodogram class #
################################################################################

setClass(
  Class = 'Periodogram',
  contains = c('Spectrum'))



################################################################################
# Object construction                                      Periodogram objects #
################################################################################

setGeneric(
  name = 'Periodogram',
  def  = function(x, ...) standardGeneric('Periodogram'))

setMethod(
  f   = 'Periodogram',
  sig = c(x = 'Waveform'),
  def = function(x) {
    # Use attributes of the Waveform x to determine Spectrum attributes.
    .nyquist     <- sampleRate(x) / 2
    .bin.width   <- sampleRate(x) / length(samples(x))
    .frequencies <- seq(from = 0, to = .nyquist, by = .bin.width)
    # Compute the periodogram ordinates.
    .ordinates   <- (1 / (N(x)*samplePeriod(x))) * abs(fft(samples(x)))^2
    .ordinates   <- .ordinates[1:length(.frequencies)]
    # Construct a new Periodogram object.
    new(Class = 'Periodogram',
        values   = .ordinates,
        binWidth = .bin.width,
        nyquist  = .nyquist)
  })

