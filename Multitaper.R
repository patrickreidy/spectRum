# Author:       Patrick Reidy
# Email:        <patrick.francis.reidy@gmail.com>
# Affiliations: The Ohio State University, Department of Linguistics
# Date:         October 22, 2014
# Purpose:      Define the Multitaper S4 class and methods for manipulating
#               objects of that class.
# Dependencies: multitaper



# Quietly load the dependencies.
suppressPackageStartupMessages(library('multitaper'))





################################################################################
# Class definition                                            Multitaper class #
################################################################################

setClass(
  Class = 'Multitaper',
  contains = c('Spectrum'),
  slots    = c(tapers = 'DPSS',
               k      = 'numeric',
               nw     = 'numeric'))



################################################################################
# Object construction                                       Multitaper objects #
################################################################################

setGeneric(
  name = 'Multitaper',
  def  = function(x, k, nw, ...) standardGeneric('Multitaper'))

setMethod(
  f   = 'Multitaper',
  sig = c(x = 'Waveform', k = 'numeric', nw = 'numeric'),
  def = function(x, k, nw) {
    # Use attributes of the Waveform x to determine Spectrum attributes.
    .nyquist     <- sampleRate(x) / 2
    .bin.width   <- sampleRate(x) / length(samples(x))
    .frequencies <- seq(from = 0, to = .nyquist, by = .bin.width)
    # Compute the DPSS tapers.
    .dpss    <- DPSS(x, k, nw)
    # Make k tapered copies of the Waveform x.
    .tapered <- apply(tapers(.dpss), 2, `*`, samples(x))
    # Compute the k eigenspectra of the tapered copies of x.
    .eigen   <- (1 / (N(x)*samplePeriod(x))) * abs(fft(.tapered))^2
    .eigen   <- .eigen[1:length(.frequencies), ]
    # Average to eigenspectra to compute the multitaper spectrum.
    .multitaper <- rowMeans(.eigen)
    # Construct a new Multitaper object.
    new(Class = 'Multitaper',
        values   = .multitaper,
        binWidth = .bin.width,
        nyquist  = .nyquist,
        k        = k,
        nw       = nw,
        tapers   = .dpss)
  })



################################################################################
# Methods                                                  Multitaper methods  #
################################################################################

#############################################################################
# @k get-method                                     MTS order parameter (k) #
#############################################################################

# k
if (! isGeneric('k'))
  setGeneric(
    name = 'k',
    def  = function(x) standardGeneric('k'))

setMethod(
  f   = 'k',
  sig = c(x = 'Multitaper'),
  def = function(x) x@k)

# nTapers
if (! isGeneric('nTapers'))
  setGeneric(
    name = 'nTapers',
    def  = function(x) standardGeneric('nTapers'))

setMethod(
  f   = 'nTapers',
  sig = c(x = 'Multitaper'),
  def = function(x) k(x))


#############################################################################
# @nw get-method                         DPSS time-bandwidth parameter (nw) #
#############################################################################

# nw
if (! isGeneric('nw'))
  setGeneric(
    name = 'nw',
    def  = function(x) standardGeneric('nw'))

setMethod(
  f   = 'nw',
  sig = c(x = 'Multitaper'),
  def = function(x) x@nw)


#############################################################################
# Peak Hz                                                           Peak Hz #
#############################################################################

# peakHz
if (! isGeneric('peakHz'))
  setGeneric(
    name = 'peakHz',
    def  = function(x, minHz, maxHz) standardGeneric('peakHz'))

setMethod(
  f   = 'peakHz',
  sig = c(x = 'Multitaper', minHz = 'numeric', maxHz = 'numeric'),
  def = function(x, minHz, maxHz) {
    .inds <- which(minHz <= frequencies(x) & frequencies(x) <= maxHz)
    .freq <- frequencies(x)[.inds]
    .vals <- values(x)[.inds]
    .freq[which.max(.vals)]
  })


#############################################################################
#  show                                                                show #
#############################################################################

setMethod(
  f   = 'show',
  sig = c(object = 'Multitaper'),
  def = function(object) {
    .spec = data.frame(Frequency = frequencies(object),
                       Amplitude = 10*log10(values(object)/max(values(object))))
    print(
      ggplot(data = .spec, aes(x = Frequency, y = Amplitude)) +
        geom_path(colour = 'black') + theme_bw() +
        xlab('Frequency (Hz)') + ylab('Amplitude (dB)')
    )
  })


#############################################################################
# @tapers get-method                                            DPSS tapers #
#############################################################################

# tapers
if (! isGeneric('tapers'))
  setGeneric(
    name = 'tapers',
    def  = function(x) standardGeneric('tapers'))

setMethod(
  f   = 'tapers',
  sig = c(x = 'Multitaper'),
  def = function(x) x@tapers)