# <p >Trace DSP, a time series processor for SAC and miniSEED data</p>

1. [Name](#)
1. [Synopsis](#synopsis)
1. [Description](#description)
1. [Options](#options)
1. [Metadata Files](#metadata-files)
1. [List Files](#list-files)
1. [Author](#author)

## <a id='synopsis'>Synopsis</a>

<pre >
tracedsp [options] file1 [file2 file3 ...]
</pre>

## <a id='description'>Description</a>

<p ><b>tracedsp</b> is a digital signal processor for time series data in SAC and miniSEED format.  For each continuous time series in the input data the specified operations will be performed in order.</p>

<p >If '-' is specified standard input will be read.  Multiple input files will be processed in the order specified.  If an input file name is prefixed with an '@' character the file is assumed to contain a list of input data files, see <i>LIST FILES</i> below.</p>

<p >When a input file is full SEED including both SEED headers and data records all of the headers will be skipped and completely unprocessed.</p>

## <a id='options'>Options</a>

<b>-V</b>

<p style="padding-left: 30px;">Print program version and exit.</p>

<b>-h</b>

<p style="padding-left: 30px;">Print program usage and exit.</p>

<b>-v</b>

<p style="padding-left: 30px;">Be more verbose.  This flag can be used multiple times ("-v -v" or "-vv") for more verbosity.</p>

<b>-s</b>

<p style="padding-left: 30px;">Print a basic summary including the number of records and the number of samples they included after processing all input records.</p>

<b>-tt </b><i>secs</i>

<p style="padding-left: 30px;">Specify a time tolerance for constructing continous trace segments. The tolerance is specified in seconds.  The default tolerance is 1/2 of the sample period.</p>

<b>-rt </b><i>diff</i>

<p style="padding-left: 30px;">Specify a sample rate tolerance for constructing continous trace segments. The tolerance is specified as the difference between two sampling rates.  The default tolerance is tested as: (abs(1-sr1/sr2) < 0.0001).</p>

<b>-ts </b><i>time</i>

<p style="padding-left: 30px;">Limit processing to miniSEED records that start after <i>time</i>. The format of the <i>time</i> arguement is: 'YYYY[,DDD,HH,MM,SS,FFFFFF]' where valid delimiters are either commas (,), colons (:) or periods (.).</p>

<b>-te </b><i>time</i>

<p style="padding-left: 30px;">Limit processing to miniSEED records that end before <i>time</i>. The format of the <i>time</i> arguement is: 'YYYY[,DDD,HH,MM,SS,FFFFFF]' where valid delimiters are either commas (,), colons (:) or periods (.).</p>

<b>-Im or -Is</b>

<p style="padding-left: 30px;">Force input data to be considered either (m) miniSEED or (s) SAC.  By default the input format is auto-detected.  This override takes effect for all files specified on the command line after the option.</p>

<b>-MSEED</b>

<p style="padding-left: 30px;">Write output as miniSEED instead of SAC.</p>

<b>-Mr </b><i>bytes</i>

<p style="padding-left: 30px;">Specify the input record length in <i>bytes</i>.  By default the length of the first miniSEED record is automatically detected, this option forces the record length.  The option is required when the input records do not contain a 1000 Blockette.</p>

<b>-Me </b><i>encoding</i>

<p style="padding-left: 30px;">Specify the data encoding format.  These encoding values are the same as those specified in the SEED 1000 Blockette.</p>

<b>-MR </b><i>bytes</i>

<p style="padding-left: 30px;">Specify the record length in <i>bytes</i> for output miniSEED data. The default record length will be 4096 bytes.</p>

<b>-ME </b><i>encoding</i>

<p style="padding-left: 30px;">Specify the <i>encoding</i> format output miniSEED data.  This can be one of:</p>
<pre style="padding-left: 30px;">
3  : 32-bit integers
4  : 32-bit floating point numbers
5  : 64-bit floating point numbers
10 : Steim-1 integer compression
11 : Stime-2 integer compression
</pre>

<b>-Sf </b><i>format</i>

<p style="padding-left: 30px;">Specify the SAC output type, the default is SAC binary in host computer byte order.  This can be one of:</p>
<pre style="padding-left: 30px;">
1 : SAC alphanumeric format
2 : SAC binary format in host computer byte order (default)
3 : SAC binary in little-endian byte order
4 : SAC binary in big-endian byte order
</pre>

<b>-m </b><i>metafile</i>

<p style="padding-left: 30px;">Specify a file containing metadata such as scaling factor, coordinates, elevation, component orientation, etc.  For each time series written any matching metadata will be added to the SAC header.  see <i>METADATA FILES</i> below.</p>

<b>-msi</b>

<p style="padding-left: 30px;">Convert any component inclination values in a metadata file from SEED (dip) to SAC convention, this is a simple matter of adding 90 degrees.</p>

<b>-c </b><i>channel</i>

<p style="padding-left: 30px;">Specify a <i>channel/component</i> for processed output data, this effectively overrides the channel component in the input data.</p>

<b>-o </b><i>outfile</i>

<p style="padding-left: 30px;">Specify the output file name.  This is optional and if not specified file names will be generated based on source parameters and data start time.  SAC files cannot contain more than a single segment of data, if this option is used and more than a single SAC segment is written as output the file will be overwritten.</p>

<b>-od </b><i>directory</i>

<p style="padding-left: 30px;">If no output file is explicitly specified the program will generate output file name(s).  This option controls which directories these file(s) will be written to.</p>

<b>-plog </b><i>logfile</i>

<p style="padding-left: 30px;">Write log entries detailing the processing applied to each segment to the specified file.</p>

<b>-CR </b><i>respfile[:start:stop]</i>

<p style="padding-left: 30px;">Convolve the time series with the response specified in the SEED RESP or StationXML file.  Optionally the starting and ending stages can be limited by specifying them as <i>start</i> and <i>stop</i>.</p>

<b>-DR </b><i>respfile[:start:stop]</i>

<p style="padding-left: 30px;">Deconvolve the response specified in the SEED RESP or StationXML file from the time series.  Optionally the starting and ending stages can be limited by specifying them as <i>start</i> and <i>stop</i>.</p>

<b>-CS </b><i>pzfile</i>

<p style="padding-left: 30px;">Convolve the time series with the response specified in the SAC poles and zeros file.</p>

<b>-DS </b><i>pzfile</i>

<p style="padding-left: 30px;">Deconvolve the response specified in the SAC poles and zeros file from the time series.</p>

<b></b><b>SPECIAL CASE</b>

<p style="padding-left: 30px;">When a deconvolution operation is directly followed by a convolution operation they will be combined into a composite transfer operation for stability and efficiency.</p>

<b>-FL </b><i>f1/f2/f3/f4</i>

<p style="padding-left: 30px;">Specify the frequency limits for deconvolution operations, frequencies are specified in Hertz.  The limits are applied as a cosine taper that scales the spectrum from 0 to 1 between f1 and f2 and from 1 to 0 between f3 and f4 in the frequency domain.</p>

<b>-FLa lowerdBdown[/upperdBdown]</b>

<p style="padding-left: 30px;">Automatically determine frequency limits for deconvolution.  A pass band is determined for all frequencies with the lower and upper corner cutoffs defined in terms of dB down from the maximum amplitude.  The upper corner definition is optional and defaults to 3 dB down.  This algorithm is designed to work with flat responses, i.e. a reponse in velocity for an instrument which is flat to velocity.  Other combinations will likely result in unsatisfactory results, i.e. a SAC PZ response in displacement for a sensor flat to velocity.</p>

<b>-W</b>

<p style="padding-left: 30px;">Prewhiten the time series data before and dewhiten after any (de)convolution operations.  Whitening is performed in the time domain.</p>

<b>-Wo </b><i>order</i>

<p style="padding-left: 30px;">Specify the prediction filter <i>order</i> used for whitening, default order is 6.</p>

<b>-Rts</b>

<p style="padding-left: 30px;">Use the total channel sensitivity when (de)convolving with responses specified in SEED RESP or StationXML.  By default evalresp calculates the total gain as the product of stage gains.  This option allows using a subset of the stages combined with the full scale sensitivity, primarily useful to calculate a response using just the poles and zeros stage with the full sensitivity (as would be done with a SAC poles and zeros file).  This is equivalent to the -ts evalresp option.</p>

<b>-Red</b>

<p style="padding-left: 30px;">Use the estimated delay in SEED RESP or StationXML files instead of the correction applied value when (de)convolving asymmetric FIR filters.  This is equivalent to the -use-estimated-delay evalresp option.</p>

<b>-Rin</b>

<p style="padding-left: 30px;">Ignore the SEED naming parameters when (de)convolving with responses specified in SEED RESP or StationXML format, i.e. use the first response found in the specified file.  By default the network, station, location and channel of the time series are matched with an appropriate response in the specified file.</p>

<b>-Ru </b><i>units</i>

<p style="padding-left: 30px;">Specify the output units for the response calculated with evalresp (de)convolution operations.  By default the output units of the response will be the input units specified in the response file.</p>

<b>-LP </b><i>frequency[/order]</i>

<p style="padding-left: 30px;">Low-pass filter the time series using an IIR filter derived from a low pass cutoff in Hertz and a filter order.  The filter \fPorder\fP can optionally be specified and defaults to 4.  The filter is applied in the forward and reverse directions to eliminate phase distortion.  The argument <b>-LP1</b> can be used to request a single pass filter, phase distortion might be present.</p>

<b>-HP </b><i>frequency[/order]</i>

<p style="padding-left: 30px;">High-pass filter the time series using an IIR filter derived from a high pass cutoff in Hertz and a filter order.  The filter \fPorder\fP can optionally be specified and defaults to 4.  The filter is applied in the forward and reverse directions to eliminate phase distortion. The argument <b>-HP1</b> can be used to request a single pass filter, phase distortion might be present.</p>

<b>-BP </b><i>frequency[/order]:frequency[/order]</i>

<p style="padding-left: 30px;">Band-pass filter the time series using an IIR filter derived from low and high pass cutoff frequencies in Hertz and filter orders.  The filter orders can optionally be specified and default to 4.  The filter is applied in the forward and reverse directions to eliminate phase distortion.  The argument <b>-BP1</b> can be used to request a single pass filter, phase distortion might be present.</p>

<b></b><b>FILTER OPERATIONS RETAIN DC OFFSET </b>

<p style="padding-left: 30px;">The <b>-LP</b>, <b>-HP</b> and <b>-BP</b> filtering operations retain the original DC offset by removing the mean value prior to filtering and restoring the mean value after the filter operation is complete.</p>

<b>-D2</b>

<p style="padding-left: 30px;">Perform a 2-point, uncentered differentiation on the time series. This results in one less sample and a time-shift of 1/2 sample period.</p>

<b>-IT</b>

<p style="padding-left: 30px;">Perform integration the time series using the trapezoidal (midpoint) method.  This results in one less sample and a time-shift of 1/2 sample period.</p>

<b>-RM</b>

<p style="padding-left: 30px;">Remove the mean from the time series.</p>

<b>-SC </b><i>factor</i>

<p style="padding-left: 30px;">Scale the time series by <i>factor</i>, i.e. multiple each data sample by <i>factor</i>.</p>

<b>-SI </b><i>factor</i>

<p style="padding-left: 30px;">Scale the time series by the inverse of <i>factor</i>, i.e. divide each data sample by <i>factor</i>.</p>

<b>-DEC </b><i>factor</i>

<p style="padding-left: 30px;">Decimate the time series by <i>factor</i> and apply an anti-alias FIR filter.  The decimation <i>factor</i> must be between 2 and 7.  The hardcoded linear-phase anti-alias filters are the same default filters used by SAC and should not disrupt the phase characteristics of the signal.</p>

<b>-TAP </b><i>width[:type]</i>

<p style="padding-left: 30px;">Apply symmetric taper of to the time series.  The taper window <i>width</i> is specified as a percent of the trace length from 0 to 0.5.  An optional window type may be specified, supported types are:</p>

<pre style="padding-left: 30px;">
HANNING (default)
HAMMING
COSINE
</pre>

<b>-POLYM </b><i>c0,c1,c2,...</i>

<p style="padding-left: 30px;">Apply a Maclaurin type polynomial to the time series using the specified coefficients as follows:</p>

<pre style="padding-left: 30px;">
output(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n
</pre>

<p >where <b>c0</b>,<b>c1</b>,<b>c2</b> .. <b>cn</b> are the coefficients and <b>x</b> is the input sample.</p>

<b>-ENV</b>

<p style="padding-left: 30px;">Calculate the envelope of the time series.  This calculation uses a Hilbert transform approximated by a time domain filter.</p>

<b>-DTRIM</b>

<p style="padding-left: 30px;">Trim each data segment to the common extents for each channel.  In other words, trim each segment to start at the latest start time of any channel and end at the earliest end time of any channel.  In different words, make all channels start and end at the same times.</p>

<b>-ROTATE E[/1],N[/2],Z[/3]:azimuth[,incidence]</b>

<p style="padding-left: 30px;">Rotate component sets.  The first three values specify the channel orientation code of the <b>east</b>, <b>north</b> and <b>vertical</b> components that should be considered a channel set.  The vertical component is optional for 2-D rotations.  Each of these codes may optionally be followed by a new code for the rotated trace.</p>

<p style="padding-left: 30px;">2-D rotations will be performed when only an <b>azimuth</b> is specified.  The <b>E</b> and <b>N</b> components will be rotated <b>azimuth</b> degrees clockwise from north.</p>

<p style="padding-left: 30px;">3-D rotations to the ray oriented LQT system will be performed when both <b>azimuth</b> and <b>incidence</b> and all three components are specified.  Assuming Z, N, and E are positive the L component will be positive along the ray defined by the azimuth and incidence, Q will be orthgonal to L in the vertical plane positive up and T will be orthogonal to both in the horizontal plane positive clockwise from the azimuth.</p>

<p style="padding-left: 30px;">For example, a simple 35 degree 2-D rotation of horizontal components:</p>

<pre style="padding-left: 30px;">
-ROTATE E,N:35
</pre>

<p >The same rotation but renaming the orientation codes to T and R:</p>

<pre >
-ROTATE E/T,N/R:35
</pre>

<p >Another example of 3-D rotation to the ray oriented LQT system 35 degrees clockwise from the original north axis and 18.8 degress of incidence with the original vertical axis:</p>

<pre >
-ROTATE E/T,N/Q,Z/L:35,18.8
</pre>

<p >If the input data is SAC and the original orientation values are set in the header they will be updated appropriately.</p>

<b>-STATS</b>

<p style="padding-left: 30px;">Calculate simple series statistics for each segment and add to the process log (see the <b>-plog</b> option), the verbose option will also cause them to be printed to stderr.  The statistics include: minimum, maximum, mean, standard deviation and RMS.</p>

## <a id='metadata-files'>Metadata Files</a>

<p >A metadata file contains a list of station parameters, some of which can be stored in SAC but not in miniSEED.  Each line in a metadata file should be a comma-separated list of parameters in the following order:</p>

<pre >
Network (KNETWK)
Station (KSTNM)
Location (KHOLE)
Channel (KCMPNM)
Latitude (STLA)
Longitude (STLO)
Elevation (STEL), in meters [not currently used by SAC]
Depth (STDP), in meters [not currently used by SAC]
Component Azimuth (CMPAZ), degrees clockwise from north
Component Incident Angle (CMPINC), degrees from vertical
Instrument Name (KINST), up to 8 characters
Scale Factor (SCALE)
Scale Frequency, unused
Scale Units, unused
Sampling rate, unused
Start time, used for matching
End time, used for matching


For example:
------------------
#net,sta,loc,chan,lat,lon,elev,depth,azimuth,SACdip,instrument,scale,scalefreq,scaleunits,samplerate,start,end
IU,ANMO,00,BH1,34.945981,-106.457133,1671,145,328,90,Geotech KS-54000,3456610000,0.02,M/S,20,2008-06-30T20:00:00,2599-12-31T23:59:59
IU,ANMO,00,BH2,34.945981,-106.457133,1671,145,58,90,Geotech KS-54000,3344370000,0.02,M/S,20,2008-06-30T20:00:00,2599-12-31T23:59:59
IU,ANMO,00,BHZ,34.945981,-106.457133,1671,145,0,0,Geotech KS-54000,3275080000,0.02,M/S,20,2008-06-30T20:00:00,2599-12-31T23:59:59
IU,ANMO,10,BH1,34.945913,-106.457122,1767.2,48.8,64,90,Guralp CMG3-T,32805600000,0.02,M/S,40,2008-06-30T20:00:00,2599-12-31T23:59:59
IU,ANMO,10,BH2,34.945913,-106.457122,1767.2,48.8,154,90,Guralp CMG3-T,32655000000,0.02,M/S,40,2008-06-30T20:00:00,2599-12-31T23:59:59
IU,ANMO,10,BHZ,34.945913,-106.457122,1767.2,48.8,0,0,Guralp CMG3-T,33067200000,0.02,M/S,40,2008-06-30T20:00:00,2599-12-31T23:59:59
------------------

As a special case '--' can be used to match a blank (space, space)
location code.
</pre>

<p >For each time series written, metadata from the first line with matching source name parameters (network, station, location and channel) and time window (if specified) will be inserted into the SAC header.  All parameters are optional except for the first four fields specifying the source name parameters.</p>

<p >Simple wildcarding: for the source name parameters that will be matched a '*' character in a field will match anything.  The BHZ metadata lines above, for example, can be (almost) summarized as:</p>

<pre >
IU,ANMO,*,BHZ,34.9459,-106.4571,1671,145,0,0,Geotech KS-54000,3456610000,0.02,M/S,20,2008-06-30T20:00:00,2599-12-31T23:59:59
</pre>

## <a id='list-files'>List Files</a>

<p >If an input file is prefixed with an '@' character the file is assumed to contain a list of file for input.  Multiple list files can be combined with multiple input files on the command line.  The last, space separated field on each line is assumed to be the file name to be read.</p>

<p >An example of a simple text list:</p>

<pre >
TA.ELFS..LHE.R.mseed
TA.ELFS..LHN.R.mseed
TA.ELFS..LHZ.R.mseed
</pre>

## <a id='author'>Author</a>

<pre >
Chad Trabant
IRIS Data Management Center

In an effort to avoid reinventing the wheel and creating new bugs many
of the core data processing routines were borrowed from other
developments including, but not limited to, SAC 2000, PQLX, Seismic
Handler and others.  Any new bugs introduced are mine.
</pre>


(man page 2018/08/21)
