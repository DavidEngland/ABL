Atmospheric Surface Flux Group Data from "Met City".

Processed data from flux group tower is given in 4 tab-delimited ASCII files.  The site, instrument, and data processing descriptions can be found in Andreas et al (1999) and Persson et al (2001).
Andreas, E.L., C. W. Fairall, P. S. Guest, and P. O. G. Persson, 1999: An overview of the SHEBA atmospheric flux program.  Preprints, Fifth Conf. on Polar Meteor. And Ocean., Jan. 10-15, 1999, Dallas, TX, 550-555.

Persson, P. Ola G., C. W. Fairall, E. Andreas, P. Guest and D. Perovich, 2001: Measurements near the               Atmospheric Surface Flux Group tower at SHEBA Part I: Site description, data processing, and accuracy estimates.  J. Geophys. Res. Submitted.
Other ASFG SHEBA publications describe additional details of specific instruments.  A list of these publications can be found on the ASFG website http://www.weather.nps.navy.mil/~psguest/sheba/

Two files contain hourly values and two files contain daily averages.  Our current version is ASFG3.0 and was completed in January 2001.  Fluxes are calculated using the observed surface pressure (at Florida) rather than an assumed constant one.  Wind direction is true wind direction accounting for the rotation of the tower during the year.  Both objective and some subjective editing has been done at various stages of the data processing.  Hourly averages were calculated as long as at least 4 10-minute periods during the hour contained 2 or more minutes of good data.  Fluxes are also eliminated when the airflow was from the ship or through the tower.  The data does include intercomparison calibrations done during the year and has been corrected based on the intercomparisons and methods described in http://www.weather.nps.navy.mil/~psguest/sheba/tower_cal/ and Persson et al 2001. 

One file (prof_file_all6_ed_hd.txt) contains data from all 5 levels on the tower as well as the radiometer and surface measurements (include the Barnes, GE and the two thermistors).  The tower data includes sensible heat and momentum fluxes as well as the RH, wind speed, temperature data.  This file also includes variances, structure functions, etc for each level.  The other file (main_file6_hd.txt) is derived from the first file by taking the median fluxes and interpolating the temperature, humidity and wind to 2.5 and 10-m.  This file does not contain any of the variances, etc, but does contain bulk estimates of latent and sensible heat fluxes at 2.5 and 10 m.  Each file has a header that explains each column of data.

There are headers in each of the files which will give you an idea of what each column is.  I will list them below so you might have an idea which file you want.

For file prof_file_all6_ed_hd.txt, the column headers are:

JD	lat	lon	Press	z1	z2	z3	z4	z5	zoph	ws1	ws2	ws3	ws4	ws5	wd1	wd2	wd3	wd4	wd5	T1	T2	T3	T4	T5	q1	q2	q3	q4	q5	rh1	rh2	rh3	rh4	rh5	rhi1	rhi2	rhi3	rhi4	rhi5	T_GE	Td_GE	T_s_epp	T_s_brns	Tsnw	Tice	Tsfc	LWd	Lwu	SWd	Swu	RR_org	RR_ncr	twr_orien	u*1	u*2	u*3	u*4	u*5	hs1	hs2	hs3	hs4	hs5	hl	ww1	ww2	ww3	ww4	ww5	sgu1	sgu2	sgu3	sgu4	sgu5	sgv1	sgv2	sgv3	sgv4	sgv5	sgw1	sgw2	sgw3	sgw4	sgw5	sgT1	sgT2	sgT3	sgT4	sgT5	cu21	cu22	cu23	cu24	cu25	cv21	cv22	cv23	cv24	cv25	cw21	cw22	cw23	cw24	cw25	cT21	cT22	cT23	cT24	cT25	No1	No2	No3	No4	No5	fl1	fl2	fl3	fl4	fl5		

with respective units:
	
UTC	deg	deg	mb	m	m	m	m	m	m	m/s	m/s	m/s	m/s	m/s	deg true	deg true	deg true	deg true	deg true	deg C	deg C	deg C	deg C	deg C	g/kg	g/kg	g/kg	g/kg	g/kg	%	%	%	%	%	%	%	%	%	%	deg C	deg C	deg C	deg C	deg C	deg C	deg C	W/m^2	W/m^2	W/m^2	W/m^2	mm/h	mm/h	deg true	m/s	m/s	m/s	m/s	m/s	W/m^2	W/m^2	W/m^2	W/m^2	W/m^2	W/m^2	m/s	m/s	m/s	m/s	m/s	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	m/s^2	deg C^2	deg C^2	deg C^2	deg C^2	deg C^2	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3	x^2/m^2/3												


For main_file6_hd.txt, the column headings are:
JJD	lat	lon	Press	ws2.5	ws10	wd2.5	wd10	T2.5	T10	q2.5	q10	rhi2.5	rhi10	T_sfc	T_s_epp	T_s_brns	Tice	LWd	Lwu	SWd	Swu	RR_org	RR_ncr	twr_orien	ww	u*	hs	hl	usb_2.5	usb_10	hsb_2.5	hsb_10	hlb_2.5	hlb_10

with the respective units:
	

UTC	deg	deg	mb	m/s	m/s	deg true	deg true	deg C	deg C	g/kg	g/kg	%	%	deg C	deg C	deg C	deg C	W/m^2	W/m^2	W/m^2	W/m^2	mm/h	mm/h	deg true	m/s	m/s	W/m^2	W/m^2	m/s	m/s	W/m^2	W/m^2	W/m^2	W/m^2

The Julian Day is based on Jan. 1, 1997 being Julian Day 1, so Jan. 1, 1998 is Julian Day 366.  The daily files contain similar parameters, and are called prof_file_davg_all6_ed_hd.txt and main_file_davg6_hd.txt.  Daily means were only calculated if at least 20 hours of data was available.  The columns "No#" and "fl#" in prof_file_all6_ed_hd.txt indicate the number of points used in the turbulence FFTs and a flag indicating whether the turbulence data passed (=0) or failed(=1) various quality control checks, respectively.  If fl=1, the wind data might be suspect, as it is also determined from the sonic anemometers.

Note that rhi is relative humidity with respect to ice recalculated from the original measurement of relative humidity with respect to water.  twr_orien is the orientation of the tower sonic boom arms, which was used to calculate the true wind direction wd.  Three surface temperature measurements are available from the General Eastern, the Eppley radiometer and the Barnes radiometer.  The Eppley is the most reliable, though there are periods when the other two are also reasonable, and one period (May) when the Eppley data may be slightly off.  Tsfc is our best estimate of the surface temperature, and is principally based on slight corrections to the Eppley temperatures and the Barnes temperatures when the Eppley was known to be wrong.  Check with one of the PIs on the latest status of the estimated errors.  The measurements of stress and sensible heat flux are the median values of the levels with "good" measurements.  Eddy correlation measurements of the latent heat flux from the Ophir instruments are included, but are quite low with respect to bulk estimates so should be used with caution.  The bulk estimates of stress, sensible and latent heat flux are calculated using a modified COARE flux algorithm that computes fluxes over the ocean or sea ice.  For ice it uses Andreas 1987 for Ch and Ce; presently it sets zo=4.5e-4 m.

Other data files are also available from the four PIs.  These data files include similar data as above but with higher temporal resolution, files of covariance and quadrature spectra for each hour and level, and sodar data.  This data is currently at various stages of development.  Two manuscripts describing this dataset and the near-surface conditions determined from this data set have been submitted for publication and can be requested from Ola Persson (ola.persson@noaa.gov).

Please let us know if you have any questions, and especially let us know if you find anything odd or an obvious error.  

Ola Persson


