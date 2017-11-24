#!/bin/csh -f

########################################################################
#
# Script: radiation_atmos_av_mon.csh
# Author: Charles Seman (reference: Steve Klein)
#
#=======================================================================
#
# ERBE obs data source:
#
#  Datasets constructed by Mark Crane:
#    /net/jjp/Datasets/obs/{lwcldfrc,swcldfrc,nswt,olr}.toa.mon.ltm.nc
#
#  have been consolidated into a compressed tar data file which is used as
#  a data source for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/ERBE/data.tar.gz
#
#=======================================================================
#
# CERES ERBE-like ES-4 obs data source:
#
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/level3_es4_table.html
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/ES4/Quality_Summaries/CER_ES4_Terra_Edition2.html
#
#  NASA CERES ERBE-like ES-4: Terra FM1 Edition2 Rev1
#  2.5 Degree Regional, Monthly (Day), Total-sky and Clear-sky
#  Averages
#
#  User Applied Revision "Rev1" applied to monthly CERES Terra
#  FM1 Edition2 variables for 5-year Climatology (3/00 - 2/05)
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/ES4/Quality_Summaries/CER_ES4_Terra_Edition2.html#user_revision
#
#  variables reference: Tables 4-4, 4-5 at
#  http://asd-www.larc.nasa.gov/ceres/collect_guide/ES4_CG.pdf
#
#  Referencing Data in Journal Articles:
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/ES4/Quality_Summaries/CER_ES4_Terra_Edition2.html#Referencing
#
#   Journal Article Reference:
#
#  Wielicki, B. A., B. R. Barkstrom, E. F. Harrison, R. B. Lee III, G. L. Smith,
#  and J. E. Cooper, 1996: Clouds and the Earth's Radiant Energy System (CERES):
#  An Earth Observing System Experiment, Bull. Amer. Meteor. Soc., 77, 853-868.
#
#   Acknowledgment:
#
#  "These data were obtained from the Atmospheric Science Data Center at
#   NASA Langley Research Center."
#
#  Acknowledgment from:
#  http://eosweb.larc.nasa.gov/GUIDE/dataset_documents/cer_es4.html#ack
#
#  "These data were obtained from the NASA Langley Research Center
#   Atmospheric Science Data Center."
#
#  Datasets:
#    CERES_ES4_rdir/lwcf.CER_ES4_Terra-FM1_Edition2.climo.nc
#    CERES_ES4_rdir/netcf.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc
#    CERES_ES4_rdir/swcf.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc
#    CERES_ES4_rdir/Total-sky/net_radiant_flux.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc
#    CERES_ES4_rdir/Total-sky/longwave_flux.CER_ES4_Terra-FM1_Edition2.climo.nc
#    CERES_ES4_rdir/Total-sky/swabs.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc
#    where CERES_ES4_rdir = /net/cjs/data/ceres/CER_ES4/Terra_Edition2/climatology/200003-200502/2.5_Degree_Regional/Monthly-Day
#
#  have been consolidated into a compressed tar data file which is used as
#  a data source for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_ES4/data.tar.gz
#
#=======================================================================
#
# CERES SRBAVG obs data source:
#
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/level3_srbavg_table.html
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/SRBAVG/Quality_Summaries/CER_SRBAVG_Edition2D_Terra_Aqua.html
#
#  NASA CERES SRBAVG Terra (CERES-FM1 or CERES-FM2) Edition2D-Rev1
#  1 Degree Monthly
#
#  User Applied Revision "Rev1" applied to monthly CERES Terra
#  (CERES-FM1 or CERES-FM2) Edition2D nonGEO and GEO SRBAVG1 variables
#  for 5-year Climatology (3/00 - 2/05)
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/SRBAVG/Quality_Summaries/CER_SRBAVG_Edition2D_Terra_Aqua.html#user_revision
#
#  Attribution - Referencing Data in Journal Articles from:
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/SRBAVG/Quality_Summaries/CER_SRBAVG_Edition2D_Terra_Aqua.html#refer
#
#   Journal Article Reference:
#
#  Wielicki, B. A., B. R. Barkstrom, E. F. Harrison, R. B. Lee III, G. L. Smith,
#  and J. E. Cooper, 1996: Clouds and the Earth's Radiant Energy System (CERES):
#  An Earth Observing System Experiment, Bull. Amer. Meteor. Soc., 77, 853-868.
#
#   Acknowledgment:
#
#  "These data were obtained from the NASA Langley Research Center EOSDIS
#   Distributed Active Archive Center."
# 
#  Acknowledgment from:
#  http://eosweb.larc.nasa.gov/GUIDE/dataset_documents/cer_srbavg.html#acknowledge
#
#  "These data were obtained from the NASA Langley Research Center
#   Atmospheric Science Data Center."
#
#  For CERES SRBAVG nonGEO, datasets:
#    /net/cjs/data/ceres/CER_SRBAVG/SRBAVG1/climatology/CER_SRBAVG1_Terra-FM1-Edition2D.nonGEO-Rev1.climo.{ctl,ieee}
#
#  have been consolidated into a compressed tar data file which is used as
#  a data source for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_SRBAVG_nonGEO/data.tar.gz
#
#  For CERES SRBAVG GEO, datasets:
#    /net/cjs/data/ceres/CER_SRBAVG/SRBAVG1/climatology/CER_SRBAVG1_Terra-FM1-Edition2D.GEO-Rev1.climo.{ctl,ieee}
#    /net/cjs/data/ceres/CER_SRBAVG/SRBAVG1/climatology/CER_SRBAVG1_Terra-XTRK-Edition2D.GEO-Rev1.climo.{ctl,ieee}
#    /net/cjs/data/ceres/CER_SRBAVG/SRBAVG1/climatology/CER_SRBAVG1_Terra-XTRK-FM1-Edition2D.GEO-Rev1.climo.{ctl,ieee}
#
#  have been consolidated into a compressed tar data file which is used as
#  a data source for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_SRBAVG_GEO/data.tar.gz
#
#=======================================================================
#
# CERES EBAF Edition1A obs data source:
#
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/level4_ebaf_table.html
#
#  NASA CERES EBAF Terra (CERES-FM1 or CERES-FM2) Edition1A
#  Energy Balanced and Filled 1 Degree Monthly (3/00 - 10/05) and
#  5-year Climatology (3/00 - 2/05)
#
#  Attribution from:
#  http://eosweb.larc.nasa.gov/PRODOCS/ceres/EBAF/Quality_Summaries/CER_EBAF_Terra_Edition1A.html
#
#   Journal Article Reference (updated using http://ams.allenpress.com/perlserv/?request=cite-builder&doi=10.1175%2F2008JCLI2637.1):
#
#  Loeb, N.G., B.A. Wielicki, D.R. Doelling, G.L. Smith, D.F. Keyes, S. Kato,
#  N. Manalo-Smith, and T. Wong, 2009: Toward Optimal Closure of the Earth's 
#  Top-of-Atmosphere Radiation Budget. J. Climate, 22, 748-766.
#
#   Acknowledgment:
#
#  "These data were obtained from the NASA Langley Research Center EOSDIS
#   Distributed Active Archive Center."
#
#  Acknowledgment from:
#  http://eosweb.larc.nasa.gov/GUIDE/dataset_documents/cer_ebaf.html#acknowledge
#
#  "These data were obtained from the NASA Langley Research Center
#   Atmospheric Science Data Center."
#
#  Datasets:
#    CERES_EBAF_dir/{clim_lwcre,clim_netcre,clim_net,clim_lwup,clim_swinc,clim_swup,clim_swcre}.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc
#    where CERES_EBAF_dir = /net/cjs/data/ceres/CERES_EBAF_TOA_Terra_Edition1A/climatology
#
#  have been consolidated into a compressed tar data file which is used as
#  a data source for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_EBAF/data.tar.gz
#
#=======================================================================
#
# CERES EBAF Edition2.6 obs data source:
#
#  CERES Product Information:
#  http://ceres.larc.nasa.gov/products.php?product=EBAF : "Browse & Subset", to:
#
#  CERES_EBAF-TOA-Terra_Ed2.6 Subsetting:
#  http://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSelection.jsp
#
#  NASA CERES EBAF Terra (CERES-FM1 or CERES-FM2) Edition2.6
#  Energy Balanced and Filled 1 Degree Monthly (3/00 - 12/10) and
#  10-year Climatology
#
#  Note, the following "Product Attribution" section was copied from the bottom of page 3 of
#  /net/cjs/data/ceres/CERES_EBAF-TOA_Terra_Ed2.6/monthly/CERES_EBAF-TOA_Terra_Ed2.6_DataProductCatalog.pdf
#  into /net/cjs/data/ceres/CERES_EBAF-TOA_Terra_Ed2.6/monthly/Product_Attribution,
#  edited there and is copied from that file here; CERES_EBAF-TOA_Terra_Ed2.6_DataProductCatalog.pdf
#  was downloaded as part of the download procedure to obtain the monthly time series dataset.
#
#  Product Attribution:
#
#  The CERES Team has gone to considerable trouble to remove major errors and to verify the quality
#  and accuracy of this data. Please specify the CERES product and version as "CERES EBAF Ed2.6"
#  and provide a reference to the following paper when you publish scientific results with the data:
#
#  Loeb, N.G., B.A. Wielicki, D.R. Doelling, G.L. Smith, D.F. Keyes, S. Kato, N. Manalo-Smith, and T.
#  Wong, 2008: Toward Optimal Closure of the Earth's Top-of-Atmosphere Radiation Budget, Journal of
#  Climate, Volume 22, Issue 3 (February 2009) pp. 748-766. doi: 10.1175/2008JCLI2637.1
#
#  Datasets:
#    CERES_EBAF_dir/monthly/CERES_EBAF-TOA_Terra_Ed2.6_Subset_200003-201012.nc
#    CERES_EBAF_dir/climate/CERES_EBAF-TOA_Terra_Ed2.6_Subset_CLIM01-CLIM12.nc
#    where CERES_EBAF_dir = /net/cjs/data/ceres/CERES_EBAF-TOA_Terra_Ed2.6
#
#  have been consolidated into a compressed tar data file which is used as
#  a data source for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_EBAF_Ed2.6/data.tar.gz
#
#=======================================================================
#
# CERES EBAF Edition2.7 obs data sources:
#
#- - - - - - - - - - - - - - - - - TOA - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF TOA Edition2.7
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 4/2013) and
#  13-year Climatology
#
#  CERES Data Products:
#  http://ceres.larc.nasa.gov/order_data.php : "EBAF-TOA", to:
#
#  CERES EBAF-TOA Product Information:
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA : "Browse & Order", to:
#
#  CERES_EBAF-TOA_Ed2.7 Subsetting and Browsing:
#  http://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSelection.jsp
#
#  the following "Product Statement" was copied from the bottom of
#  http://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                             Product Statement
#                          - CERES_EBAF-TOA_Ed2.7 -
#
#  Product Description:
#
#  The CERES_EBAF-TOA_Ed2.7 provides monthly and climatological averages of
#  clear-sky fluxes, all-sky fluxes, and cloud radiative effect (CRE) fluxes at
#  TOA, that are energy balanced to the ocean heat storage term and clear-sky
#  spatially filled. More information can be obtained here: CERES-EBAF-TOA
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA.
#
#  Product Data Quality Summary:
#
#  The EBAF-TOA data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF_Ed2.7_DQS.pdf
#  provides more information of the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-TOA Ed2.7" and provide a reference to the following
#  paper when you publish scientific results with the data:
#
#    Loeb, N.G., B.A. Wielicki, D.R. Doelling, G.L. Smith, D.F. Keyes, S. Kato,
#    N. Manalo-Smith, and T. Wong, 2009: Toward Optimal Closure of the Earth's
#    Top-of-Atmosphere Radiation Budget. Journal of Climate, Volume 22, Issue 3
#    (February 2009) pp. 748-766. doi: 10.1175/2008JCLI2637.1
#
#- - - - - - - - - - - - - - - - - SFC - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF Surface Edition2.7
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 9/2012) and
#  12-year Climatology
#
#  CERES Data Products:
#  http://ceres.larc.nasa.gov/order_data.php : "EBAF-Surface", to:
#
#  CERES EBAF-Surface Product Information:
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-Surface: "Browse & Order", to:
#
#  CERES_EBAF-Surface_Ed2.7 Subsetting and Browsing:
#  http://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSFCSelection.jsp
#
#  the following "Product Statement" was copied from the bottom of
#  http://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                             Product Statement
#                        - CERES_EBAF-Surface_Ed2.7 -
#
#  Product Description:
#
#  The CERES_EBAF-Surface_Ed2.7 product provides monthly and climatological
#  averages of computed clear-sky fluxes, all-sky fluxes, and cloud radiative
#  effect (CRE) fluxes at surface, consistent with the CERES EBAF-TOA fluxes. More
#  information can be obtained here: CERES-EBAF-Surface
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-Surface.
#
#  Product Data Quality Summary:
#
#  The EBAF-Surface data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF-Surface_Ed2.7_DQS.pdf
#  provides more information of the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-Surface Ed2.7" and provide a reference to the
#  following paper when you publish scientific results with the data:
#
#    Kato, S., N. G. Loeb, F. G. Rose, D. R. Doelling, D. A. Rutan, T. E.
#    Caldwell, L. Yu, and R. A. Weller, 2013: Surface irradiances consistent with
#    CERES-derived top-of-atmosphere shortwave and longwave irradiances. Journal
#    of Climate, Volume 26, 2719-2740. doi: 10.1175/JCLI-D-12-00436.1
#
#  Datasets:
#    CERES_dir/CERES_EBAF-Surface_Ed2.7/monthly/CERES_EBAF-Surface_Ed2.7_Subset_200003-201209.nc
#    CERES_dir/CERES_EBAF-Surface_Ed2.7/climate/CERES_EBAF-Surface_Ed2.7_Subset_CLIM01-CLIM12.nc
#    CERES_dir/CERES_EBAF-TOA_Ed2.7/monthly/CERES_EBAF-TOA_Ed2.7_Subset_200003-201304.nc
#    CERES_dir/CERES_EBAF-TOA_Ed2.7/climate/CERES_EBAF-TOA_Ed2.7_Subset_CLIM01-CLIM12.nc
#    where CERES_dir = /net2/cjs/data/ceres
#
#  have been consolidated into a compressed tar data file which is used as
#  a data source for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_EBAF_Ed2.7/data.tar.gz
#
#=======================================================================
#
# CERES EBAF Edition2.8 v1 obs data sources:
#
#- - - - - - - - - - - - - - - - - TOA - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF TOA Edition2.8
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 3/2015) and
#  15-year Climatology
#
#  CERES Data Products:
#  http://ceres.larc.nasa.gov/order_data.php : "EBAF-TOA", to:
#
#  CERES EBAF-TOA Product Information:
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA : "Browse & Order", to:
#
#  CERES_EBAF-TOA_Ed2.8 Subsetting and Browsing:
#  http://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSelection.jsp
#
#  the following "Product Statement" is an edited version of a copy
#  taken from the bottom of http://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                             Product Statement
#                          - CERES_EBAF-TOA_Ed2.8 -
#
#  Product Description:
#
#  The CERES_EBAF-TOA_Ed2.8 provides monthly and climatological averages of
#  clear-sky fluxes, all-sky fluxes, and cloud radiative effect (CRE) fluxes at
#  TOA, that are energy balanced to the ocean heat storage term and clear-sky
#  spatially filled. More information can be obtained here: CERES-EBAF-TOA
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA .
#
#  Product Data Quality Summary:
#
#  The EBAF-TOA data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF_Ed2.8_DQS.pdf
#  provides more information of the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-TOA Ed2.8" and provide a reference to the following
#  paper when you publish scientific results with the data:
#
#    Loeb, N.G., B.A. Wielicki, D.R. Doelling, G.L. Smith, D.F. Keyes, S. Kato,
#    N. Manalo-Smith, and T. Wong, 2009: Toward Optimal Closure of the Earth's
#    Top-of-Atmosphere Radiation Budget. Journal of Climate, Volume 22, Issue 3
#    (February 2009) pp. 748-766. doi: 10.1175/2008JCLI2637.1
#
#- - - - - - - - - - - - - - - - - SFC - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF Surface Edition2.8
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 2/2015) and
#  15-year Climatology
#
#  CERES Data Products:
#  http://ceres.larc.nasa.gov/order_data.php : "EBAF-Surface", to:
#
#  CERES EBAF-Surface Product Information:
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-Surface: "Browse & Order", to:
#
#  CERES_EBAF-Surface_Ed2.8 Subsetting and Browsing:
#  http://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSFCSelection.jsp
#
#  the following "Product Statement" is an edited version of a copy
#  taken from the bottom of http://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                             Product Statement
#                        - CERES_EBAF-Surface_Ed2.8 -
#
#  Product Description:
#
#  The CERES_EBAF-Surface_Ed2.8 product provides monthly and climatological
#  averages of computed clear-sky fluxes, all-sky fluxes, and cloud radiative
#  effect (CRE) fluxes at surface, consistent with the CERES EBAF-TOA fluxes. More
#  information can be obtained here: CERES-EBAF-Surface
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-Surface .
#
#  Product Data Quality Summary:
#
#  The EBAF-Surface data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF-Surface_Ed2.8_DQS.pdf
#  provides more information of the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-Surface Ed2.8" and provide a reference to the
#  following paper when you publish scientific results with the data:
#
#    Kato, S., N. G. Loeb, F. G. Rose, D. R. Doelling, D. A. Rutan, T. E.
#    Caldwell, L. Yu, and R. A. Weller, 2013: Surface irradiances consistent with
#    CERES-derived top-of-atmosphere shortwave and longwave irradiances. Journal
#    of Climate, Volume 26, 2719-2740. doi: 10.1175/JCLI-D-12-00436.1
#
#  Datasets:
#    CERES_dir/CERES_EBAF-Surface_Ed2.8/v1/monthly/CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc
#    CERES_dir/CERES_EBAF-Surface_Ed2.8/v1/climate/CERES_EBAF-Surface_Ed2.8_Subset_CLIM01-CLIM12.nc
#    CERES_dir/CERES_EBAF-TOA_Ed2.8/v1/monthly/CERES_EBAF-TOA_Ed2.8_Subset_200003-201503.nc
#    CERES_dir/CERES_EBAF-TOA_Ed2.8/v1/climate/CERES_EBAF-TOA_Ed2.8_Subset_CLIM01-CLIM12.nc
#    where CERES_dir = /net2/cjs/data/ceres
#
#  have been consolidated into a compressed tar data file for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_EBAF_Ed2.8/v1/data.tar.gz
#
#=======================================================================
#
# CERES EBAF Edition2.8 v2 obs data sources:
#
#- - - - - - - - - - - - - - - - - TOA - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF TOA Edition2.8
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 7/2016) and
#  16-year Climatology
#
#  CERES Data Products:
#  http://ceres.larc.nasa.gov/order_data.php : "EBAF-TOA", to:
#
#  CERES EBAF-TOA Product Information:
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA : "Browse & Subset", to:
#
#  CERES_EBAF-TOA_Ed2.8 Subsetting and Browsing:
#  http://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSelection.jsp
#
#  the following "Product Statement" is an edited version of a copy
#  taken from the bottom of https://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                             Product Statement
#                          - CERES_EBAF-TOA_Ed2.8 -
#
#  Product Description:
#
#  The CERES_EBAF-TOA_Ed2.8 provides monthly and climatological averages of
#  clear-sky fluxes, all-sky fluxes, and cloud radiative effect (CRE) fluxes at
#  TOA, that are energy balanced to the ocean heat storage term and clear-sky
#  spatially filled. More information can be obtained here: CERES-EBAF-TOA
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA .
#
#  Product Data Quality Summary:
#
#  The EBAF-TOA data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF_Ed2.8_DQS.pdf
#  provides more information of the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-TOA Ed2.8" and provide a reference to the following
#  paper when you publish scientific results with the data:
#
#      Loeb, N.G., B.A. Wielicki, D.R. Doelling, G.L. Smith, D.F. Keyes, S. Kato,
#      N. Manalo-Smith, and T. Wong, 2009: Toward Optimal Closure of the Earth's
#      Top-of-Atmosphere Radiation Budget. Journal of Climate, Volume 22, Issue 3
#      (February 2009) pp. 748-766. doi: 10.1175/2008JCLI2637.1
#
#- - - - - - - - - - - - - - - - - SFC - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF Surface Edition2.8
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 2/2016) and
#  16-year Climatology
#
#  CERES Data Products:
#  http://ceres.larc.nasa.gov/order_data.php : "EBAF-Surface", to:
#
#  CERES EBAF-Surface Product Information:
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-Surface: "Browse & Subset", to:
#
#  CERES_EBAF-Surface_Ed2.8 Subsetting and Browsing:
#  http://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSFCSelection.jsp
#
#  the following "Product Statement" is an edited version of a copy
#  taken from the bottom of https://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                             Product Statement
#                        - CERES_EBAF-Surface_Ed2.8 -

#  Product Description:
#
#  The CERES_EBAF-Surface_Ed2.8 product provides monthly and climatological
#  averages of computed clear-sky fluxes, all-sky fluxes, and cloud radiative
#  effect (CRE) fluxes at surface, consistent with the CERES EBAF-TOA fluxes. More
#  information can be obtained here: CERES-EBAF-Surface
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-Surface .
#
#  Product Data Quality Summary:
#
#  The EBAF-Surface data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF-Surface_Ed2.8_DQS.pdf
#  provides more information of the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-Surface Ed2.8" and provide a reference to the
#  following paper when you publish scientific results with the data:
#
#      Kato, S., N. G. Loeb, F. G. Rose, D. R. Doelling, D. A. Rutan, T. E.
#      Caldwell, L. Yu, and R. A. Weller, 2013: Surface irradiances consistent with
#      CERES-derived top-of-atmosphere shortwave and longwave irradiances. Journal
#      of Climate, Volume 26, 2719-2740. doi: 10.1175/JCLI-D-12-00436.1
#
#  Datasets:
#    CERES_dir/CERES_EBAF-Surface_Ed2.8/v2/monthly/CERES_EBAF-Surface_Ed2.8_Subset_200003-201602.nc
#    CERES_dir/CERES_EBAF-Surface_Ed2.8/v2/climate/CERES_EBAF-Surface_Ed2.8_Subset_CLIM01-CLIM12.nc
#    CERES_dir/CERES_EBAF-TOA_Ed2.8/v2/monthly/CERES_EBAF-TOA_Ed2.8_Subset_200003-201607.nc
#    CERES_dir/CERES_EBAF-TOA_Ed2.8/v2/climate/CERES_EBAF-TOA_Ed2.8_Subset_CLIM01-CLIM12.nc
#    where CERES_dir = /net2/cjs/data/ceres
#
#  have been consolidated into a compressed tar data file for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_EBAF_Ed2.8/v2/data.tar.gz
#
#=======================================================================
#
# CERES EBAF Edition4.0 v1 obs data sources:
#
#- - - - - - - - - - - - - - - - - TOA - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF TOA Edition4.0
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 1/2017) and
#  7/2005 to 6/2015 Climatology
#
#  CERES Data Products:
#  https://ceres.larc.nasa.gov/order_data.php : "EBAF-TOA", to:
#
#  CERES EBAF-TOA Product Information:
#  https://ceres.larc.nasa.gov/products.php?product=EBAF-TOA : "Browse & Subset", to:
#
#  CERES EBAF-TOA Ordering Page:
#  https://ceres.larc.nasa.gov/products-info.php?product=EBAF-TOA : "Browse & Subset", to:
#
#  CERES Subsetting and Browsing (CERES_EBAF-TOA_Ed4.0 Subsetting and Browsing):
#  https://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAF4Selection.jsp
#
#  the following "Product Statement" is an edited version of a copy
#  taken from the bottom of https://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                                 Product Statement
#                              - CERES_EBAF-TOA_Ed4.0 -
#
#  Product Description:
#
#  The CERES_EBAF-TOA_Ed4.0 provides monthly and climatological averages of
#  clear-sky fluxes, all-sky fluxes and clouds, and cloud radiative effect (CRE)
#  fluxes at TOA, that are energy balanced to the ocean heat storage term and
#  clear-sky spatially filled. More information can be obtained here: CERES-EBAF-TOA
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA .
#
#  Product Data Quality Summary:
#
#  The EBAF-TOA data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF_Ed4.0_DQS.pdf
#  provides more information on the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-TOA Ed4.0" and provide a reference to the following
#  paper when you publish scientific results with the data:
#
#      Loeb, N.G., B.A. Wielicki, D.R. Doelling, G.L. Smith, D.F. Keyes, S. Kato,
#      N. Manalo-Smith, and T. Wong, 2009: Toward Optimal Closure of the Earth's
#      Top-of-Atmosphere Radiation Budget. Journal of Climate, Volume 22, Issue 3
#      (February 2009) pp. 748-766. doi: 10.1175/2008JCLI2637.1
#
#- - - - - - - - - - - - - - - - - SFC - - - - - - - - - - - - - - - - -
#
#  NASA CERES EBAF Surface Edition4.0
#  Energy Balanced and Filled 1 Degree Monthly (3/2000 - 2/2016) and
#  7/2005 to 6/2015 Climatology
#
#  CERES Data Products:
#  https://ceres.larc.nasa.gov/order_data.php : "EBAF-Surface", to:
#
#  CERES EBAF-Surface Product Information:
#  https://ceres.larc.nasa.gov/products.php?product=EBAF-Surface : "Browse & Subset", to:
#
#  CERES EBAF-Surface Ordering Page:
#  https://ceres.larc.nasa.gov/products-info.php?product=EBAF-Surface : "Browse & Subset", to:
#
#  CERES_EBAF-Surface_Ed4.0 Subsetting and Browsing:
#  https://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAFSFC4Selection.jsp
#
#  the following "Product Statement" is an edited version of a copy
#  taken from the bottom of https://ceres-tool.larc.nasa.gov/ord-tool/srbavg
#  during the "Get Data" part of the data order process:
#
#                                 Product Statement
#                           - CERES_EBAF-Surface_Ed4.0 -
#
#  Product Description:
#
#  The CERES_EBAF-Surface_Ed4.0 product provides monthly and climatological
#  averages of computed clear-sky fluxes, all-sky fluxes, and cloud radiative
#  effect (CRE) fluxes at surface, consistent with the CERES EBAF-TOA fluxes. More
#  information can be obtained here: CERES-EBAF-Surface
#  http://ceres.larc.nasa.gov/products.php?product=EBAF-Surface .
#
#  Product Data Quality Summary:
#
#  The EBAF-Surface data quality summary
#  http://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF-Surface_Ed4.0_DQS.pdf
#  provides more information of the content.
#
#  Product Attribution:
#
#  The CERES Team has made considerable efforts to remove major errors and to
#  verify the quality and accuracy of this data. Please specify the CERES product
#  and version as "CERES EBAF-Surface Ed4.0" and provide a reference to the
#  following paper when you publish scientific results with the data:
#
#      Kato, S., N. G. Loeb, F. G. Rose, D. R. Doelling, D. A. Rutan, T. E.
#      Caldwell, L. Yu, and R. A. Weller, 2013: Surface irradiances consistent with
#      CERES-derived top-of-atmosphere shortwave and longwave irradiances. Journal
#      of Climate, Volume 26, 2719-2740. doi: 10.1175/JCLI-D-12-00436.1
#
#  Datasets:
#    CERES_dir/CERES_EBAF-Surface_Ed4.0/v1/monthly/CERES_EBAF-Surface_Ed4.0_Subset_200003-201602.nc
#    CERES_dir/CERES_EBAF-Surface_Ed4.0/v1/climate/CERES_EBAF-Surface_Ed4.0_Subset_CLIM01-CLIM12.nc
#    CERES_dir/CERES_EBAF-TOA_Ed4.0/v1/monthly/CERES_EBAF-TOA_Ed4.0_Subset_200003-201701.nc
#    CERES_dir/CERES_EBAF-TOA_Ed4.0/v1/climate/CERES_EBAF-TOA_Ed4.0_Subset_CLIM01-CLIM12.nc
#    where CERES_dir = /net2/cjs/data/ceres
#
#  have been consolidated into a compressed tar data file for this processor:
#
#  $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/CERES_EBAF_Ed4.0/v1/data.tar.gz
#
#=======================================================================
#
# Model data source ("atmos" for GFDL named variables):
#
#   pp/atmos/av/monthly_Xyr
#
#=======================================================================
#
# Output:
#
#   $out_dir/atmos_${yr1}_${yr2}
#
#=======================================================================
#
# Sample frepp usage (http://www.gfdl.noaa.gov/fre#Post-Processing):
#
# <component type="atmos" zInterp="era40" start="1981" source="atmos_month" cubicToLatLon="90,144">
#    <timeAverage source="monthly" interval="20yr">
#       <analysis script="$FRE_ANALYSIS_HOME/cjs/stub/radiation_atmos_av_mon.csh"/>
#    </timeAverage>
# </component>
#
########################################################################

echo ; uname -a ; echo

echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo $0 $*
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

# make sure the required modules are loaded from parent shell
if (`gfdl_platform` == "hpcs-csc") then
   source $MODULESHOME/init/csh
   echo
   module list
   echo
   module use -a /home/fms/local/modulefiles
   set error = 0
   if (`echo '$M="fre";@M=grep/^$M\/.*/,split/:/,$ENV{"LOADEDMODULES"};if($M[0]=~/^$M\/(.*)$/){print $1}'|perl` == "") then
      echo "ERROR: module fre not loaded."; @ error++
   endif
   if (`echo '$M="fre-analysis";@M=grep/^$M\/.*/,split/:/,$ENV{"LOADEDMODULES"};if($M[0]=~/^$M\/(.*)$/){print $1}'|perl` == "") then
      echo "ERROR: module fre-analysis not loaded."; @ error++
   endif
   if ($error > 0) exit 1
else
   echo "Running on `gfdl_platform` to run script on data prepared for publishing"
endif

# Parse argument list
set argv = (`getopt Ui:o:d:y: $*`)

while ("$argv[1]" != "--")
    switch ($argv[1])
        case -U:
            set NoCompress; breaksw
        case -i:
            set in_data_dir = $argv[2]; shift argv; breaksw
        case -o:
            set out_dir = $argv[2]; shift argv; breaksw
        case -d:
            set descriptor = $argv[2]; shift argv; breaksw
        case -y:
            set years = $argv[2]; shift argv; breaksw
    endsw
    shift argv
end
shift argv

# argument error checking

if ( $#argv == 1 || $#argv == 12 ) then
   set in_data_file = ( $argv )
else
   echo "ERROR: number of input files must be 1 or 12."
   set help
endif

if ( "$in_data_dir" == "" ) then
   echo "ERROR: no argument given for input directory."
   set help
else
   if ( ! -d "$in_data_dir" ) then
      echo "ERROR: argument given for input directory is not a directory: "$in_data_dir
      set help
   endif
endif

if ( "$out_dir" == "" ) then
   echo "ERROR: no argument given for output directory."
   set help
endif

if ( "$descriptor" == "" ) then
   echo "ERROR: no argument given for descriptor."
   set help
endif

if ( $?years ) then
   set yrs = `echo $years | sed -e "s/,/ /"`
   if ( $#yrs != 2 ) then
      echo "ERROR: invalid entry for years."
      set help
   else
      set yr1 = $yrs[1]; set yr2 = $yrs[2]
   endif
else
   echo "ERROR: no entry for years."
   set help
endif

if ( $?help ) then
   echo
   echo "USAGE:  $0:t -i idir -o odir -d desc -y yrs files...."
   echo
   exit 1
endif

# use yrs format for figures directory - if supplied

if ( "$yr1" != "" && "$yr2" != "" ) then
   set out_dir = ${out_dir}/atmos_${yr1}_${yr2}
else
   echo 'yr1 yr2: '$yr1 $yr2
   exit 1
endif

echo 'in_data_dir: '$in_data_dir
echo 'in_data_file: '$in_data_file
echo 'out_dir: '$out_dir
echo 'descriptor: '$descriptor
echo 'yr1 yr2: '$yr1 $yr2
if ( $?NoCompress ) then
  echo 'Requested no compression for postscript files...'
else
  echo 'Requested compression for postscript files...'
endif

########################################################################

echo
set WORK_DIR = $TMPDIR
if ( ! -d $WORK_DIR ) mkdir -p $WORK_DIR || exit 1
cd $WORK_DIR || exit 1
echo 'current working directory: '$cwd
if ( $cwd == $TMPDIR ) then
  echo 'remove any old files in '$cwd
  wipetmp
endif
echo

echo '-----------------------------------------------------------------'
echo '                         model data'
echo '-----------------------------------------------------------------'

echo
echo $descriptor'...'
echo

#  model file work directory...

set mod_file_work_dir = $WORK_DIR/INPUT/$descriptor
mkdir -p $mod_file_work_dir || exit 1

cd $mod_file_work_dir || exit 1

echo 'working in directory:'
echo $cwd

echo
echo 'dmget and copy model file(s) to current work directory:'
echo
set filelist = ( )
foreach file ( $in_data_file )
  set file_path = $in_data_dir/$file
  set filelist = ( $filelist $file_path )
  ls $file_path
end
dmget $filelist || exit 1
echo
gcp $filelist . || exit 1
echo
echo 'finished dmget and file copy to current work directory:'
echo $cwd
echo

# for each model input file, check and for each variable that contains
# a "missing_value" attribute but no "_FillValue" attribute, rename the
# variable's "missing_value" attribute to an "_FillValue" attribute

foreach file ( $in_data_file )
  set missing_value_vars = ( `ncdump -h $file | grep missing_value | cut -d: -f1` )
  set _FillValue_vars = ( `ncdump -h $file | grep _FillValue | cut -d: -f1` )
  foreach var ( $missing_value_vars )
    set grep_result = ( `ncdump -h $file | grep "${var}:" | grep _FillValue` )
    if ( $#grep_result == 0 ) then	# no _FillValue; rename variable's missing_value to _FillValue
      echo "ncrename -a ${var}@missing_value,_FillValue ${file}"
      ncrename -a ${var}@missing_value,_FillValue $file || exit 1
    endif
    unset grep_result
  end
  if ( $#missing_value_vars > 0 && $#missing_value_vars > $#_FillValue_vars ) echo
end

set varlist = (olr,olr_clr,swdn_toa,swup_toa,swup_toa_clr,lwdn_sfc,lwdn_sfc_clr,lwup_sfc,lwup_sfc_clr,swdn_sfc,swdn_sfc_clr,swup_sfc,swup_sfc_clr)

# Bruce Wyman's technique to identify if there are CMIP named variables

if (`ncdump -h $in_data_file[1] | perl -e '$d=join"",<stdin>;if($d=~/\t\w+ rlut\(.+\)/){print "1"}else{print "0"}'`) then
  #
  # execute CMIP named variables
  #
  echo 'attempting to add radiation variables "'$varlist'" into'
  echo $in_data_file
  echo 'using the CMIP named variables'
  foreach file ( $in_data_file )
    echo '...'$file
    foreach var ( {$varlist} )
      if ( -e tmp.nc ) rm -f tmp.nc
      switch ( $var )
        case olr:
          ncdump -h $file | grep 'rlut(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rlut $file -o tmp.nc
            ncrename -v rlut,olr tmp.nc
            ncks -A -v olr tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "olr" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case olr_clr:
          ncdump -h $file | grep 'rlutcs(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rlutcs $file -o tmp.nc
            ncrename -v rlutcs,olr_clr tmp.nc
            ncks -A -v olr_clr tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "olr_clr" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case swdn_toa:
          ncdump -h $file | grep 'rsdt(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rsdt $file -o tmp.nc
            ncrename -v rsdt,swdn_toa tmp.nc
            ncks -A -v swdn_toa tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "swdn_toa" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case swup_toa:
          ncdump -h $file | grep 'rsut(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rsut $file -o tmp.nc
            ncrename -v rsut,swup_toa tmp.nc
            ncks -A -v swup_toa tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "swup_toa" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case swup_toa_clr:
          ncdump -h $file | grep 'rsutcs(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rsutcs $file -o tmp.nc
            ncrename -v rsutcs,swup_toa_clr tmp.nc
            ncks -A -v swup_toa_clr tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "swup_toa_clr" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case lwdn_sfc:
          ncdump -h $file | grep 'rlds(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rlds $file -o tmp.nc
            ncrename -v rlds,lwdn_sfc tmp.nc
            ncks -A -v lwdn_sfc tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "lwdn_sfc" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case lwdn_sfc_clr:
          ncdump -h $file | grep 'rldscs(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rldscs $file -o tmp.nc
            ncrename -v rldscs,lwdn_sfc_clr tmp.nc
            ncks -A -v lwdn_sfc_clr tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "lwdn_sfc_clr" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case lwup_sfc:
          ncdump -h $file | grep 'rlus(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rlus $file -o tmp.nc
            ncrename -v rlus,lwup_sfc tmp.nc
            ncks -A -v lwup_sfc tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "lwup_sfc" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case lwup_sfc_clr:
          #
          # note: "rluscs" does not exist for atmos_cmip; according to Bruce Wyman, use "rlus":
          # "LW up at the surface is the same for full sky and clear sky."
          #
          ncdump -h $file | grep 'rlus(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rlus $file -o tmp.nc
            ncrename -v rlus,lwup_sfc_clr tmp.nc
            ncks -A -v lwup_sfc_clr tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "lwup_sfc_clr" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case swdn_sfc:
          ncdump -h $file | grep 'rsds(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rsds $file -o tmp.nc
            ncrename -v rsds,swdn_sfc tmp.nc
            ncks -A -v swdn_sfc tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "swdn_sfc" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case swdn_sfc_clr:
          ncdump -h $file | grep 'rsdscs(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rsdscs $file -o tmp.nc
            ncrename -v rsdscs,swdn_sfc_clr tmp.nc
            ncks -A -v swdn_sfc_clr tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "swdn_sfc_clr" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case swup_sfc:
          ncdump -h $file | grep 'rsus(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rsus $file -o tmp.nc
            ncrename -v rsus,swup_sfc tmp.nc
            ncks -A -v swup_sfc tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "swup_sfc" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
        case swup_sfc_clr:
          ncdump -h $file | grep 'rsuscs(' >& /dev/null ; set var_status = $status
          if ( $var_status == 0 ) then
            ncks -a -c -v rsuscs $file -o tmp.nc
            ncrename -v rsuscs,swup_sfc_clr tmp.nc
            ncks -A -v swup_sfc_clr tmp.nc -o $file
            if ( $status != 0 ) then
              echo 'problem adding variable "swup_sfc_clr" into: '$file
              echo
              exit 1
            endif
          endif
          breaksw
      endsw
    end
  end
  echo 'finished.'
  echo
  set varlist = (${varlist},rlut,rlutcs,rsdt,rsut,rsutcs,rlds,rldscs,rlus,rsds,rsdscs,rsus,rsuscs)
endif

set climo_file = climatology.nc

echo 'extract "'$varlist'" from:'
echo $in_data_file
echo 'for a monthly climatology input data file "'$climo_file'"'

if ( $#in_data_file == 1 ) then
  # single monthly climo file (assumed to have 12 records)
  ncks -v $varlist -a -c -F $in_data_file[1] $climo_file || exit 1
else
  # 12 monthly climo files (cat together)
  ncrcat -v $varlist $in_data_file $climo_file || exit 1
endif

if ( ! -e $climo_file ) then
  echo "ERROR: file $climo_file does not exist."
  exit 1
endif

set mod_climo_file = $mod_file_work_dir/$climo_file

echo
echo 'local model radiation climatology file:'
echo
ls -l $mod_climo_file || exit 1

########################################################################

#  local output root directory

set local_out_rdir = $WORK_DIR/OUTPUT
mkdir -p $local_out_rdir || exit 1

#  do work in local RUN directory

set local_run_dir = $WORK_DIR/RUN
mkdir -p $local_run_dir || exit 1

cd $local_run_dir || exit 1

unset do_obs_dmget

if ( `echo "$FRE_ANALYSIS_ARCHIVE" | grep -i archive | wc -l` == 1 ) then
  set do_obs_dmget
endif

echo
echo '-----------------------------------------------------------------'
echo '                 make plots using GrADS...'
echo '-----------------------------------------------------------------'
echo
echo 'working in directory:'
echo $cwd
echo

#  get GrADS scripts...

cp $FRE_ANALYSIS_HOME/cjs/shared/gs/*.gs . || exit 1

set do_stats_file = true

if ( $do_stats_file == true ) then
  cp $FRE_ANALYSIS_HOME/cjs/shared/gs/do_stats_file/{plotmeanmod,plotmeanobs,plotsd}.gs . || exit 1
  set stats_file = $local_out_rdir/stats.txt
else
  set stats_file = 'none'
endif

#  get plot file scripts which are used to make GrADS plot files for
#  radiation variables...

cp $FRE_ANALYSIS_HOME/cjs/code/radiation_atmos_av_mon/plot_file_scripts/*_plot_file.csh . || exit 1

echo 'make XDF file for "'$mod_climo_file'" calendar for GrADS v2:'

#  (NOLEAP calendar not available in GrADS v2)...

set mod_xdf_file = ${climo_file:r}.xdf
if ( -e $mod_xdf_file ) rm -f $mod_xdf_file || exit 1

cp $FRE_ANALYSIS_HOME/cjs/shared/bin/{make_GrADS_xdf,grads_starting_date,calculate_days_in_year} . || exit 1

echo
echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
$cwd/make_GrADS_xdf $mod_climo_file $mod_xdf_file || exit 1
echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
echo

echo 'determine the model grid parameters for "'$mod_climo_file'":'

cp $FRE_ANALYSIS_HOME/cjs/shared/bin/get_model_horizontal_grid_info.csh . || exit 1

echo
echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
$cwd/get_model_horizontal_grid_info.csh $mod_climo_file || exit 1
echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
echo

source model_horizontal_grid_info

echo
echo 'files in current work directory:'
echo $cwd
echo
ls -lt * || exit 1
echo

#  plot annotation for model years titles...

set title_mod_years = '('$yr1'-'$yr2')'

#  make plots for radiation variables for each observed radiation source...
#
#For the AM4 documentation paper, the following is used.
set james = 'on'
if($james == 'on') then
set obs_list = "CERES_EBAF_TOA_Ed2.8"
else
set obs_list = "ERBE CERES_ES4 CERES_SRBAVG_nonGEO CERES_SRBAVG_GEO CERES_EBAF CERES_EBAF_Ed2.6 CERES_EBAF_TOA_Ed2.7 CERES_EBAF_SFC_Ed2.7 CERES_EBAF_ATM_Ed2.7 CERES_EBAF_TOA_Ed2.8 CERES_EBAF_SFC_Ed2.8 CERES_EBAF_ATM_Ed2.8 CERES_EBAF_TOA_Ed4.0 CERES_EBAF_SFC_Ed4.0 CERES_EBAF_ATM_Ed4.0"
endif 

foreach obs_source ($obs_list)
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo '                    obs: '$obs_source
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

#  obs data source

switch ( $obs_source )
  case CERES_EBAF_TOA_Ed2.7:
  case CERES_EBAF_SFC_Ed2.7:
  case CERES_EBAF_ATM_Ed2.7:
    set obs_source_sdir = CERES_EBAF_Ed2.7
    set ceres_ebaf_toa_ed2p7_nc = CERES_EBAF-TOA_Ed2.7_Subset_200003-201304.nc
    set ceres_ebaf_sfc_ed2p7_nc = CERES_EBAF-Surface_Ed2.7_Subset_200003-201209.nc
    breaksw
  case CERES_EBAF_TOA_Ed2.8:
  case CERES_EBAF_SFC_Ed2.8:
  case CERES_EBAF_ATM_Ed2.8:
    set obs_source_sdir = CERES_EBAF_Ed2.8/v2
    set ceres_ebaf_toa_ed2p8_nc = CERES_EBAF-TOA_Ed2.8_Subset_200003-201607.nc
    set ceres_ebaf_sfc_ed2p8_nc = CERES_EBAF-Surface_Ed2.8_Subset_200003-201602.nc
    breaksw
  case CERES_EBAF_TOA_Ed4.0:
  case CERES_EBAF_SFC_Ed4.0:
  case CERES_EBAF_ATM_Ed4.0:
    set obs_source_sdir = CERES_EBAF_Ed4.0/v1
    set ceres_ebaf_toa_ed4p0_nc = CERES_EBAF-TOA_Ed4.0_Subset_200003-201701.nc
    set ceres_ebaf_sfc_ed4p0_nc = CERES_EBAF-Surface_Ed4.0_Subset_200003-201602.nc
    breaksw
  default:
    set obs_source_sdir = $obs_source
    breaksw
endsw

#  ...source data

# move the follwing to frepp template script, so the user can set it at a high level 
set obs_data_file = $obs_data_file/$obs_source_sdir/data.tar.gz # setting it as environment variable in the template script. $FRE_ANALYSIS_ARCHIVE/cjs/radiation_atmos_av_mon/$obs_source_sdir/data.tar.gz

#  ...local storage directory

set obs_file_work_dir = $WORK_DIR/INPUT/$obs_source_sdir
mkdir -p $obs_file_work_dir || exit 1

echo
echo 'obs compressed tar file:'
echo $obs_data_file
echo
if ( $?do_obs_dmget ) then
  echo 'dmget obs compressed tar file...'
  dmget $obs_data_file || exit 1
  echo '...finished'
  echo
endif
echo 'copy obs compressed tar file to directory:'
echo $obs_file_work_dir
gcp $obs_data_file $obs_file_work_dir || exit 1
echo '...finished'

#  plot files local storage directory

set local_plot_file_dir = $local_run_dir/plot_file/$obs_source
mkdir -p $local_plot_file_dir || exit 1

set obs_file_project = $local_plot_file_dir/obs_file_project

switch ( $obs_source )
  case ERBE:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract ERBE obs files from compressed tar file...'
    echo
    set obs_files = ( {lwcldfrc,swcldfrc,nswt,olr}.toa.mon.ltm.nc )
    tar -xzf data.tar.gz $obs_files || exit 1
    foreach file ( $obs_files )
      if ( ! -e $file ) then
        echo '...'$file' does not exist, exiting script'
        exit 1
      endif
    end
    echo 'ERBE obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_files || exit 1
    echo
    echo 'set ERBE_lwcldfrc_nc = lwcldfrc.toa.mon.ltm.nc' >! $obs_file_project
    echo 'set ERBE_swcldfrc_nc = swcldfrc.toa.mon.ltm.nc' >> $obs_file_project
    echo 'set ERBE_nswt_nc     =     nswt.toa.mon.ltm.nc' >> $obs_file_project
    echo 'set ERBE_olr_nc      =      olr.toa.mon.ltm.nc' >> $obs_file_project
    #
    #  obs identification for plot file name
    #
    set obs_plot_file = erbe
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 144  obs_dlon = 2.5  obs_lon1 =   0.  obs_lon2 = 357.5
    set obs_nlat =  73  obs_dlat = 2.5  obs_lat1 = -90.  obs_lat2 =  90.
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = none  obs_tf = none
    #
    #  plot annotation
    #
    set obs_title_src = 'ERBE'
    set obs_title_yrs = '(2/85-1/89)'
    set obs_title_dsc = ''
    set obs_title_att = ''
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcf swcf netcf olr swabs netrad )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa/erbe
    breaksw
  case CERES_ES4:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES ES4 obs files from compressed tar file...'
    echo
    set obs_files = ( lwcf.CER_ES4_Terra-FM1_Edition2.climo.nc \
                      netcf.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc \
                      swcf.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc \
                      Total-sky/net_radiant_flux.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc \
                      Total-sky/longwave_flux.CER_ES4_Terra-FM1_Edition2.climo.nc \
                      Total-sky/swabs.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc )
    tar -xzf data.tar.gz $obs_files || exit 1
    foreach file ( $obs_files )
      if ( ! -e $file ) then
        echo '...'$file' does not exist, exiting script'
        exit 1
      endif
    end
    echo 'CERES ES4 obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_files || exit 1
    echo
    echo 'set CERES_ES4_lwcf_nc   = lwcf.CER_ES4_Terra-FM1_Edition2.climo.nc'                            >! $obs_file_project
    echo 'set CERES_ES4_netcf_nc  = netcf.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc'                      >> $obs_file_project
    echo 'set CERES_ES4_swcf_nc   = swcf.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc'                       >> $obs_file_project
    echo 'set CERES_ES4_netrad_nc = Total-sky/net_radiant_flux.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc' >> $obs_file_project
    echo 'set CERES_ES4_lw_nc     = Total-sky/longwave_flux.CER_ES4_Terra-FM1_Edition2.climo.nc'         >> $obs_file_project
    echo 'set CERES_ES4_swabs_nc  = Total-sky/swabs.CER_ES4_Terra-FM1_Edition2_Rev1.climo.nc'            >> $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_es4
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 144  obs_dlon = 2.5  obs_lon1 =   1.25  obs_lon2 = 358.75
    set obs_nlat =  72  obs_dlat = 2.5  obs_lat1 = -88.75  obs_lat2 =  88.75
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = none  obs_tf = none
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES ERBE-like ES-4'
    set obs_title_yrs = '(3/00-2/05)'
    set obs_title_dsc = 'CERES ERBE-like ES-4: Terra "'varInst'" "'varEditionRev'"'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_ES4/Terra-FM1_Edition2/Product_Statements.txt'
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcf swcf netcf olr swabs netrad )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa/ceres_old_ed
    breaksw
  case CERES_SRBAVG_nonGEO:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES SRBAVG nonGEO obs files from compressed tar file...'
    echo
    set obs_files = "CER_SRBAVG1_Terra-FM1-Edition2D.nonGEO-Rev1.climo.{ctl,ieee}"
    tar -xzf data.tar.gz $obs_files || exit 1
    foreach file ( $obs_files )
      if ( ! -e $file ) then
        echo '...'$file' does not exist, exiting script'
        exit 1
      endif
    end
    echo
    echo 'convert CERES SRBAVG nonGEO ctl file for local use...'
    echo
    set ctl_file = CER_SRBAVG1_Terra-FM1-Edition2D.nonGEO-Rev1.climo.ctl
    set ieee_file = CER_SRBAVG1_Terra-FM1-Edition2D.nonGEO-Rev1.climo.ieee
    grep $ieee_file $ctl_file
    if ( $status != 0 ) then
      echo 'did not find reference to '$ieee_file' in '$ctl_file
      exit 1
    endif
    set old_ctl_file = old_CER_SRBAVG1_Terra-FM1-Edition2D.nonGEO-Rev1.climo.ctl
    mv $ctl_file $old_ctl_file
    echo "dset $obs_file_work_dir/$ieee_file" > $ctl_file
    grep -v dset $old_ctl_file >> $ctl_file
    echo
    echo 'CERES SRBAVG nonGEO obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_files || exit 1
    echo
    echo 'set CERES_SRBAVG_nonGEO_ctl = CER_SRBAVG1_Terra-FM1-Edition2D.nonGEO-Rev1.climo.ctl' >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_srbavg_nongeo
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = none  obs_tf = none
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES SRBAVG nonGEO'
    set obs_title_yrs = '(3/00-2/05)'
    set obs_title_dsc = 'CERES SRBAVG nonGEO: Terra "'varInst'" SRBAVG1 "'varEditionRev'"'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_SRBAVG/Product_Statements.txt'
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcf swcf netcf olr swabs netrad )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa/ceres_old_ed
    breaksw
  case CERES_SRBAVG_GEO:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES SRBAVG GEO obs files from compressed tar file...'
    echo
    set obs_files = ( CER_SRBAVG1_Terra-FM1-Edition2D.GEO-Rev1.climo.{ctl,ieee} \
                      CER_SRBAVG1_Terra-XTRK-Edition2D.GEO-Rev1.climo.{ctl,ieee} \
                      CER_SRBAVG1_Terra-XTRK-FM1-Edition2D.GEO-Rev1.climo.{ctl,ieee} )
    tar -xzf data.tar.gz $obs_files || exit 1
    foreach file ( $obs_files )
      if ( ! -e $file ) then
        echo '...'$file' does not exist, exiting script'
        exit 1
      endif
    end
    echo
    echo 'convert CERES SRBAVG GEO ctl files for local use...'
    echo
    foreach ctl_file ( *.ctl )
      set ieee_file = $ctl_file:r.ieee
      grep $ieee_file $ctl_file
      if ( $status != 0 ) then
        echo 'did not find reference to '$ieee_file' in '$ctl_file
        exit 1
      endif
      set old_ctl_file = old_$ctl_file
      mv $ctl_file $old_ctl_file
      echo "dset $obs_file_work_dir/$ieee_file" > $ctl_file
      grep -v dset $old_ctl_file >> $ctl_file
    end
    echo
    echo 'CERES SRBAVG GEO obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_files || exit 1
    echo
    echo 'set CERES_SRBAVG_GEO_FM1_ctl      = CER_SRBAVG1_Terra-FM1-Edition2D.GEO-Rev1.climo.ctl'      >! $obs_file_project
    echo 'set CERES_SRBAVG_GEO_XTRK_ctl     = CER_SRBAVG1_Terra-XTRK-Edition2D.GEO-Rev1.climo.ctl'     >> $obs_file_project
    echo 'set CERES_SRBAVG_GEO_XTRK_FM1_ctl = CER_SRBAVG1_Terra-XTRK-FM1-Edition2D.GEO-Rev1.climo.ctl' >> $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_srbavg_geo
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = none  obs_tf = none
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES SRBAVG GEO'
    set obs_title_yrs = '(3/00-2/05)'
    set obs_title_dsc = 'CERES SRBAVG GEO: Terra "'varInst'" SRBAVG1 "'varEditionRev'"'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_SRBAVG/Product_Statements.txt'
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcf swcf netcf olr swabs netrad )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa/ceres_old_ed
    breaksw
  case CERES_EBAF:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF Edition1A obs files from compressed tar file...'
    echo
    set obs_var = clim_lwcre,clim_netcre,clim_net,clim_lwup,clim_swinc,clim_swup,clim_swcre
    set obs_files = ( {$obs_var}.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc )
    tar -xzf data.tar.gz $obs_files || exit 1
    foreach file ( $obs_files )
      if ( ! -e $file ) then
        echo '...'$file' does not exist, exiting script'
        exit 1
      endif
    end
    echo 'CERES EBAF Edition1A obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_files || exit 1
    echo
    echo 'set CERES_EBAF_lwcre_nc  =  clim_lwcre.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc' >! $obs_file_project
    echo 'set CERES_EBAF_netcre_nc = clim_netcre.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc' >> $obs_file_project
    echo 'set CERES_EBAF_net_nc    =    clim_net.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc' >> $obs_file_project
    echo 'set CERES_EBAF_lwup_nc   =   clim_lwup.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc' >> $obs_file_project
    echo 'set CERES_EBAF_swinc_nc  =  clim_swinc.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc' >> $obs_file_project
    echo 'set CERES_EBAF_swup_nc   =   clim_swup.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc' >> $obs_file_project
    echo 'set CERES_EBAF_swcre_nc  =  clim_swcre.CERES_EBAF_TOA_Terra_Edition1A_200003-200510.nc' >> $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = none  obs_tf = none
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF'
    set obs_title_yrs = '(3/00-2/05)'
    set obs_title_dsc = 'CERES EBAF Terra Edition1A'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_EBAF_TOA_Terra_Edition1A/Product_Statements.txt'
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre olr swabs netrad )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa/ceres_old_ed
    breaksw
  case CERES_EBAF_Ed2.6:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF Edition2.6 obs monthly time series file from compressed tar file...'
    echo
    set obs_file = CERES_EBAF-TOA_Terra_Ed2.6_Subset_200003-201012.nc
    tar -xzf data.tar.gz $obs_file || exit 1
    if ( ! -e $obs_file ) then
      echo '...'$obs_file' does not exist, exiting script'
      exit 1
    endif
    echo
    echo 'CERES EBAF Edition2.6 obs file in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file || exit 1
    echo
    echo "set CERES_EBAF_Ed2p6_nc = $obs_file" >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_ed2.6
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 120
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF Ed2.6'
    set obs_title_yrs = '(3/00-2/10)'
    set obs_title_dsc = 'CERES EBAF Terra Edition2.6'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_EBAF-TOA_Terra_Ed2.6/Product_Statements.txt'
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre olr swabs netrad )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir  = $local_out_rdir/Seman.radiation_toa/ceres_old_ed
    set local_rad_out_rdir2 = $local_out_rdir/Klein.radiation
    breaksw
  case CERES_EBAF_TOA_Ed2.7:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF TOA Edition2.7 obs monthly time series file from compressed tar file...'
    echo
    set obs_file = $ceres_ebaf_toa_ed2p7_nc
    tar -xzf data.tar.gz $obs_file || exit 1
    if ( ! -e $obs_file ) then
      echo '...'$obs_file' does not exist, exiting script'
      exit 1
    endif
    echo
    echo 'CERES EBAF TOA Edition2.7 obs file in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file || exit 1
    echo
    echo "set CERES_EBAF_TOA_Ed2p7_nc = $obs_file" >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_toa_ed2.7
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 156
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF TOA Ed2.7'
    set obs_title_yrs = '(3/00-2/13)'
    set obs_title_dsc = 'CERES EBAF TOA Edition2.7'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_EBAF_Ed2.7/Product_Statements.txt'
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre olr swabs netrad olr_clr swabs_clr )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa/ceres_old_ed
    breaksw
  case CERES_EBAF_SFC_Ed2.7:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF Surface Edition2.7 obs monthly time series file from compressed tar file...'
    echo
    set obs_file = $ceres_ebaf_sfc_ed2p7_nc
    tar -xzf data.tar.gz $obs_file || exit 1
    if ( ! -e $obs_file ) then
      echo '...'$obs_file' does not exist, exiting script'
      exit 1
    endif
    echo
    echo 'CERES EBAF Surface Edition2.7 obs file in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file || exit 1
    echo
    echo "set CERES_EBAF_SFC_Ed2p7_nc = $obs_file" >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_sfc_ed2.7
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 144
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF Surface Ed2.7'
    set obs_title_yrs = '(3/00-2/12)'
    set obs_title_dsc = 'CERES EBAF Surface Edition2.7'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_EBAF_Ed2.7/Product_Statements.txt'
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre netlw netsw netrad \
                         lwdn lwdn_clr lwup lwup_clr \
                         swdn swdn_clr swup swup_clr )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_sfc/ceres_old_ed
    breaksw
  case CERES_EBAF_ATM_Ed2.7:
    cd $obs_file_work_dir || exit 1
    set obs_file_toa = $ceres_ebaf_toa_ed2p7_nc
    if ( ! -e $obs_file_toa ) then
      echo
      echo 'extract CERES EBAF TOA Edition2.7 obs monthly time series file from compressed tar file...'
      echo
      tar -xzf data.tar.gz $obs_file_toa || exit 1
      if ( ! -e $obs_file_toa ) then
        echo '...'$obs_file_toa' does not exist, exiting script'
        exit 1
      endif
    endif
    set obs_file_sfc = $ceres_ebaf_sfc_ed2p7_nc
    if ( ! -e $obs_file_sfc ) then
      echo
      echo 'extract CERES EBAF Surface Edition2.7 obs monthly time series file from compressed tar file...'
      echo
      tar -xzf data.tar.gz $obs_file_sfc || exit 1
      if ( ! -e $obs_file_sfc ) then
        echo '...'$obs_file_sfc' does not exist, exiting script'
        exit 1
      endif
    endif
    echo
    echo 'CERES EBAF Edition2.7 obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file_toa $obs_file_sfc || exit 1
    echo
    echo "set CERES_EBAF_TOA_Ed2p7_nc = $obs_file_toa" >! $obs_file_project
    echo "set CERES_EBAF_SFC_Ed2p7_nc = $obs_file_sfc" >> $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_atm_ed2.7
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 144
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF TOA & Surface Ed2.7'
    set obs_title_yrs = '(3/00-2/12)'
    set obs_title_dsc = 'CERES EBAF TOA & Surface Edition2.7'
    set obs_title_att = '/net/cjs/data/ceres/attribution/CERES_EBAF_Ed2.7/Product_Statements.txt'
    #
    #  diagnostic variables (available: netlw_diff, lwup_diff, netsw_diff, swup_diff, swcre_ratio)
    #
    set diag_varlist = ( netlw_diff netsw_diff swcre_ratio )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_atm/ceres_old_ed
    breaksw
  case CERES_EBAF_TOA_Ed2.8:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF TOA Edition2.8 obs monthly time series file from compressed tar file...'
    echo
    set obs_file = $ceres_ebaf_toa_ed2p8_nc
    tar -xzf data.tar.gz $obs_file || exit 1
    if ( ! -e $obs_file ) then
      echo '...'$obs_file' does not exist, exiting script'
      exit 1
    endif
    echo
    echo 'CERES EBAF TOA Edition2.8 obs file in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file || exit 1
    echo
    echo "set CERES_EBAF_TOA_Ed2p8_nc = $obs_file" >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_toa_ed2.8
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 180
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF TOA Ed2.8'
    set obs_title_yrs = '(3/00-2/15)'
    set obs_title_dsc = 'CERES EBAF TOA Edition2.8'
    set obs_title_att = "/net/cjs/data/ceres/attribution/${obs_source_sdir}/Product_Statements.txt"
    #
    # set parameters for a 16-year climatology if environmental variable is set
    #
    if ( $?do_CERES_EBAF_Ed2p8_16yr ) then
      set obs_dts = 0  obs_tf = 192
      set obs_title_yrs = '(3/00-2/16)'
    endif
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre olr swabs netrad olr_clr swabs_clr )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa
    breaksw
  case CERES_EBAF_SFC_Ed2.8:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF Surface Edition2.8 obs monthly time series file from compressed tar file...'
    echo
    set obs_file = $ceres_ebaf_sfc_ed2p8_nc
    tar -xzf data.tar.gz $obs_file || exit 1
    if ( ! -e $obs_file ) then
      echo '...'$obs_file' does not exist, exiting script'
      exit 1
    endif
    echo
    echo 'CERES EBAF Surface Edition2.8 obs file in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file || exit 1
    echo
    echo "set CERES_EBAF_SFC_Ed2p8_nc = $obs_file" >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_sfc_ed2.8
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 180
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF Surface Ed2.8'
    set obs_title_yrs = '(3/00-2/15)'
    set obs_title_dsc = 'CERES EBAF Surface Edition2.8'
    set obs_title_att = "/net/cjs/data/ceres/attribution/${obs_source_sdir}/Product_Statements.txt"
    #
    # set parameters for a 16-year climatology if environmental variable is set
    #
    if ( $?do_CERES_EBAF_Ed2p8_16yr ) then
      set obs_dts = 0  obs_tf = 192
      set obs_title_yrs = '(3/00-2/16)'
    endif
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre netlw netsw netrad \
                         lwdn lwdn_clr lwup lwup_clr \
                         swdn swdn_clr swup swup_clr )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_sfc
    breaksw
  case CERES_EBAF_ATM_Ed2.8:
    cd $obs_file_work_dir || exit 1
    set obs_file_toa = $ceres_ebaf_toa_ed2p8_nc
    if ( ! -e $obs_file_toa ) then
      echo
      echo 'extract CERES EBAF TOA Edition2.8 obs monthly time series file from compressed tar file...'
      echo
      tar -xzf data.tar.gz $obs_file_toa || exit 1
      if ( ! -e $obs_file_toa ) then
        echo '...'$obs_file_toa' does not exist, exiting script'
        exit 1
      endif
    endif
    set obs_file_sfc = $ceres_ebaf_sfc_ed2p8_nc
    if ( ! -e $obs_file_sfc ) then
      echo
      echo 'extract CERES EBAF Surface Edition2.8 obs monthly time series file from compressed tar file...'
      echo
      tar -xzf data.tar.gz $obs_file_sfc || exit 1
      if ( ! -e $obs_file_sfc ) then
        echo '...'$obs_file_sfc' does not exist, exiting script'
        exit 1
      endif
    endif
    echo
    echo 'CERES EBAF Edition2.8 obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file_toa $obs_file_sfc || exit 1
    echo
    echo "set CERES_EBAF_TOA_Ed2p8_nc = $obs_file_toa" >! $obs_file_project
    echo "set CERES_EBAF_SFC_Ed2p8_nc = $obs_file_sfc" >> $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_atm_ed2.8
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 180
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF TOA & Surface Ed2.8'
    set obs_title_yrs = '(3/00-2/15)'
    set obs_title_dsc = 'CERES EBAF TOA & Surface Edition2.8'
    set obs_title_att = "/net/cjs/data/ceres/attribution/${obs_source_sdir}/Product_Statements.txt"
    #
    # set parameters for a 16-year climatology if environmental variable is set
    #
    if ( $?do_CERES_EBAF_Ed2p8_16yr ) then
      set obs_dts = 0  obs_tf = 192
      set obs_title_yrs = '(3/00-2/16)'
    endif
    #
    #  diagnostic variables (available: netlw_diff, lwup_diff, netsw_diff, swup_diff, swcre_ratio)
    #
    set diag_varlist = ( netlw_diff netsw_diff swcre_ratio )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_atm
    breaksw
  case CERES_EBAF_TOA_Ed4.0:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF TOA Edition4.0 obs monthly time series file from compressed tar file...'
    echo
    set obs_file = $ceres_ebaf_toa_ed4p0_nc
    tar -xzf data.tar.gz $obs_file || exit 1
    if ( ! -e $obs_file ) then
      echo '...'$obs_file' does not exist, exiting script'
      exit 1
    endif
    echo
    echo 'CERES EBAF TOA Edition4.0 obs file in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file || exit 1
    echo
    echo "set CERES_EBAF_TOA_Ed4p0_nc = $obs_file" >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_toa_ed4.0
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 180
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF TOA Ed4.0'
    set obs_title_yrs = '(3/00-2/15)'
    set obs_title_dsc = 'CERES EBAF TOA Edition4.0'
    set obs_title_att = "/net/cjs/data/ceres/attribution/${obs_source_sdir}/Product_Statements.txt"
    #
    # set parameters for a 16-year climatology if environmental variable is set
    #
    if ( $?do_CERES_EBAF_Ed4p0_16yr ) then
      set obs_dts = 0  obs_tf = 192
      set obs_title_yrs = '(3/00-2/16)'
    endif
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre olr swabs netrad olr_clr swabs_clr )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_toa
    breaksw
  case CERES_EBAF_SFC_Ed4.0:
    cd $obs_file_work_dir || exit 1
    echo
    echo 'extract CERES EBAF Surface Edition4.0 obs monthly time series file from compressed tar file...'
    echo
    set obs_file = $ceres_ebaf_sfc_ed4p0_nc
    tar -xzf data.tar.gz $obs_file || exit 1
    if ( ! -e $obs_file ) then
      echo '...'$obs_file' does not exist, exiting script'
      exit 1
    endif
    echo
    echo 'CERES EBAF Surface Edition4.0 obs file in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file || exit 1
    echo
    echo "set CERES_EBAF_SFC_Ed4p0_nc = $obs_file" >! $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_sfc_ed4.0
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 180
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF Surface Ed4.0'
    set obs_title_yrs = '(3/00-2/15)'
    set obs_title_dsc = 'CERES EBAF Surface Edition4.0'
    set obs_title_att = "/net/cjs/data/ceres/attribution/${obs_source_sdir}/Product_Statements.txt"
    #
    # set parameters for a 16-year climatology if environmental variable is set
    #
    if ( $?do_CERES_EBAF_Ed4p0_16yr ) then
      set obs_dts = 0  obs_tf = 192
      set obs_title_yrs = '(3/00-2/16)'
    endif
    #
    #  diagnostic variables
    #
    set diag_varlist = ( lwcre swcre netcre netlw netsw netrad \
                         lwdn lwdn_clr lwup lwup_clr \
                         swdn swdn_clr swup swup_clr )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_sfc
    breaksw
  case CERES_EBAF_ATM_Ed4.0:
    cd $obs_file_work_dir || exit 1
    set obs_file_toa = $ceres_ebaf_toa_ed4p0_nc
    if ( ! -e $obs_file_toa ) then
      echo
      echo 'extract CERES EBAF TOA Edition4.0 obs monthly time series file from compressed tar file...'
      echo
      tar -xzf data.tar.gz $obs_file_toa || exit 1
      if ( ! -e $obs_file_toa ) then
        echo '...'$obs_file_toa' does not exist, exiting script'
        exit 1
      endif
    endif
    set obs_file_sfc = $ceres_ebaf_sfc_ed4p0_nc
    if ( ! -e $obs_file_sfc ) then
      echo
      echo 'extract CERES EBAF Surface Edition4.0 obs monthly time series file from compressed tar file...'
      echo
      tar -xzf data.tar.gz $obs_file_sfc || exit 1
      if ( ! -e $obs_file_sfc ) then
        echo '...'$obs_file_sfc' does not exist, exiting script'
        exit 1
      endif
    endif
    echo
    echo 'CERES EBAF Edition4.0 obs files in current work directory:'
    echo
    echo $cwd
    ls -lt $obs_file_toa $obs_file_sfc || exit 1
    echo
    echo "set CERES_EBAF_TOA_Ed4p0_nc = $obs_file_toa" >! $obs_file_project
    echo "set CERES_EBAF_SFC_Ed4p0_nc = $obs_file_sfc" >> $obs_file_project
    #
    #  obs identification for plot file
    #
    set obs_plot_file = ceres_ebaf_atm_ed4.0
    #
    #  obs grid anchor points for regridding model data to obs grid
    #
    set obs_nlon = 360  obs_dlon = 1.0  obs_lon1 =   0.5  obs_lon2 = 359.5
    set obs_nlat = 180  obs_dlat = 1.0  obs_lat1 = -89.5  obs_lat2 =  89.5
    #
    #  starting offset time increment and final time level for obs climatology
    #
    set obs_dts = 0  obs_tf = 180
    #
    #  plot annotation
    #
    set obs_title_src = 'CERES EBAF TOA & Surface Ed4.0'
    set obs_title_yrs = '(3/00-2/15)'
    set obs_title_dsc = 'CERES EBAF TOA & Surface Edition4.0'
    set obs_title_att = "/net/cjs/data/ceres/attribution/${obs_source_sdir}/Product_Statements.txt"
    #
    # set parameters for a 16-year climatology if environmental variable is set
    #
    if ( $?do_CERES_EBAF_Ed4p0_16yr ) then
      set obs_dts = 0  obs_tf = 192
      set obs_title_yrs = '(3/00-2/16)'
    endif
    #
    #  diagnostic variables (available: netlw_diff, lwup_diff, netsw_diff, swup_diff, swcre_ratio)
    #
    set diag_varlist = ( netlw_diff netsw_diff swcre_ratio )
    #
    #  local radiation output root directory
    #
    set local_rad_out_rdir = $local_out_rdir/Seman.radiation_atm
    breaksw
endsw

#  model plot title

set mod_title = $descriptor
set mod_title_yrs = "$title_mod_years"

#  determine if regridding of model, obs, or both will be done using model
#  and obs grid parameters; set regrid flags and define plot annotation

#  the goal is to regrid to a coarser resolution grid; this could be either
#  model, obs, or a grid based on a nominal model grid resolution (nominal
#  model grid case would be chosen if the regridding would be to the model
#  grid but the model grid is non-uniform; in this case both model and obs
#  would be regrid to the nominal model grid)

if ( $mod_lon1 == $obs_lon1 && $mod_lat1 == $obs_lat1 && \
     $mod_dlon == $obs_dlon && $mod_dlat == $obs_dlat ) then
  set regrid_dlon = none
  set regrid_lon1 = none
  set regrid_dlat = none
  set regrid_lat1 = none
  set mod_grid_type = mod ; set mod_grid_dsc = '(model grid)'
  set obs_grid_type = obs ; set obs_grid_dsc = '(obs grid)'
else
  if ( $mod_dlon == regrid_to_mod && $mod_dlat == regrid_to_mod ) then
    set regrid_dlon = $mod_dlon
    set regrid_lon1 = $mod_lon1
    set regrid_dlat = $mod_dlat
    set regrid_lat1 = $mod_lat1
    set mod_grid_type = mod           ; set mod_grid_dsc = '(model grid)'
    set obs_grid_type = regrid_to_mod ; set obs_grid_dsc = '(regrid to model grid)'
  else if ( $mod_dlon == regrid_to_obs && $mod_dlat == regrid_to_obs ) then
    set regrid_dlon = $obs_dlon
    set regrid_lon1 = $obs_lon1
    set regrid_dlat = $obs_dlat
    set regrid_lat1 = $obs_lat1
    set mod_grid_type = regrid_to_obs ; set mod_grid_dsc = '(regrid to obs grid)'
    set obs_grid_type = obs           ; set obs_grid_dsc = '(obs grid)'
  else if ( $mod_dlon == non_uniform_grid || $mod_dlat == non_uniform_grid ) then
    #
    #  define a nominal model grid (with uniform grid spacing) for further checks
    #
    if ( $mod_dlon == non_uniform_grid ) then
      set mod_regrid_dlon = `echo "360 / $mod_nlon" | bc -l` || exit 1
      set mod_regrid_lon1 = `echo "$mod_regrid_dlon / 2" | bc -l` || exit 1
    else
      set mod_regrid_dlon = $mod_dlon
      set mod_regrid_lon1 = $mod_lon1
    endif
    if ( $mod_dlat == non_uniform_grid ) then
      set mod_regrid_dlat = `echo "180 / $mod_nlat" | bc -l` || exit 1
      set mod_regrid_lat1 = `echo "-90 + $mod_regrid_dlat / 2" | bc -l` || exit 1
    else
      set mod_regrid_dlat = $mod_dlat
      set mod_regrid_lat1 = $mod_lat1
    endif
    #
    #  use nominal model grid (with uniform grid spacing)
    #  and obs grid spacings to determine regridding
    #
    set dlon_diff = `echo "$mod_regrid_dlon - $obs_dlon" | bc -l` || exit 1
    set dlat_diff = `echo "$mod_regrid_dlat - $obs_dlat" | bc -l` || exit 1
    if ( `echo "($dlon_diff + $dlat_diff) > 0" | bc -l` == 1 ) then
      set regrid_dlon = $mod_regrid_dlon
      set regrid_lon1 = $mod_regrid_lon1
      set regrid_dlat = $mod_regrid_dlat
      set regrid_lat1 = $mod_regrid_lat1
      set mod_grid_type = regrid_to_nominal_mod ; set mod_grid_dsc = '(regrid to nominal mod grid)'
      set obs_grid_type = regrid_to_nominal_mod ; set obs_grid_dsc = '(regrid to nominal mod grid)'
    else
      set regrid_dlon = $obs_dlon
      set regrid_lon1 = $obs_lon1
      set regrid_dlat = $obs_dlat
      set regrid_lat1 = $obs_lat1
      set mod_grid_type = regrid_to_obs ; set mod_grid_dsc = '(regrid to obs grid)'
      set obs_grid_type = obs           ; set obs_grid_dsc = '(obs grid)'
    endif
  else
    #
    #  use model and obs grid spacings to determine regridding
    #
    set dlon_diff = `echo "$mod_dlon - $obs_dlon" | bc -l` || exit 1
    set dlat_diff = `echo "$mod_dlat - $obs_dlat" | bc -l` || exit 1
    if ( `echo "($dlon_diff + $dlat_diff) > 0" | bc -l` == 1 ) then
      set regrid_dlon = $mod_dlon
      set regrid_lon1 = $mod_lon1
      set regrid_dlat = $mod_dlat
      set regrid_lat1 = $mod_lat1
      set mod_grid_type = mod           ; set mod_grid_dsc = '(model grid)'
      set obs_grid_type = regrid_to_mod ; set obs_grid_dsc = '(regrid to model grid)'
    else
      set regrid_dlon = $obs_dlon
      set regrid_lon1 = $obs_lon1
      set regrid_dlat = $obs_dlat
      set regrid_lat1 = $obs_lat1
      set mod_grid_type = regrid_to_obs ; set mod_grid_dsc = '(regrid to obs grid)'
      set obs_grid_type = obs           ; set obs_grid_dsc = '(obs grid)'
    endif
  endif
endif

cd $local_run_dir || exit 1

#  make local GrADS "plot" file for each variable

foreach var ( $diag_varlist )

set var_plot_file = $local_plot_file_dir/$var.${obs_plot_file}.plot
if ( -e $var_plot_file ) rm -f $var_plot_file || exit 1

set var_plot_file_project = $local_plot_file_dir/$var.plot_file_project
if ( -e $var_plot_file_project ) rm -f $var_plot_file_project || exit 1

cat << EOF > $var_plot_file_project
set diag_var = $var
set plot_file = $var_plot_file
set obs_source = $obs_source
set obs_in_data_dir = $obs_file_work_dir
source $obs_file_project
set obs_nlon = $obs_nlon
set obs_nlat = $obs_nlat
set obs_dts = $obs_dts
set obs_tf = $obs_tf
set mod_file = $mod_climo_file
set mod_xdf_file = $mod_xdf_file
set mod_nlon = $mod_nlon
set mod_nlat = $mod_nlat
set obs_grid_type = $obs_grid_type
set mod_grid_type = $mod_grid_type
set regrid_dlon = $regrid_dlon
set regrid_lon1 = $regrid_lon1
set regrid_dlat = $regrid_dlat
set regrid_lat1 = $regrid_lat1
set GrADS_version = GrADS_v2
EOF

switch ( $var )
  case lwcf:
  case lwcre:
    $cwd/lwcf_plot_file.csh $var_plot_file_project || exit 1
    breaksw
  case swcf:
  case swcre:
  case swcf_ratio:
  case swcre_ratio:
    $cwd/swcf_plot_file.csh $var_plot_file_project || exit 1
    breaksw
  case netcf:
  case netcre:
    $cwd/netcf_plot_file.csh $var_plot_file_project || exit 1
    breaksw
  case olr:
  case olr_clr:
  case netlw:
  case netlw_diff:
  case lwdn:
  case lwdn_clr:
  case lwup:
  case lwup_clr:
  case lwup_diff:
    $cwd/lwrad_plot_file.csh $var_plot_file_project || exit 1
    breaksw
  case swabs:
  case swabs_clr:
  case netsw:
  case netsw_diff:
  case swdn:
  case swdn_clr:
  case swup:
  case swup_clr:
  case swup_diff:
    $cwd/swrad_plot_file.csh $var_plot_file_project || exit 1
    breaksw
  case netrad:
    $cwd/netrad_plot_file.csh $var_plot_file_project || exit 1
    breaksw
endsw

end # var

cd $local_run_dir || exit 1

#  make GrADS plotting script

set two_back_slashes = '\\'
set three_back_slashes = '\\\'

if ( -e plot_3panel_radiation.gs ) rm -f plot_3panel_radiation.gs || exit 1
cat << EOF > plot_3panel_radiation.gs
function panel3(args)

* 1,2 prefix of model seasonal field variable name and its area average
* 3,4 prefix of obs seasonal field variable name and its area average
* 5 the variable name for the purpose of set levels
* 6 - 9 are the lon/lat to calculate statistics
* 10 is the variable unit appearing on the figures' title (if specified)
* 11 is the variable instrument (if specified)
* 12 is the variable edition revision (if specified)

mprefix    = subwrd(args,1)
maveprefix = subwrd(args,2)
oprefix    = subwrd(args,3)
oaveprefix = subwrd(args,4)
varname    = subwrd(args,5)
_lon1 = subwrd(args,6)
_lon2 = subwrd(args,7)
_lat1 = subwrd(args,8)
_lat2 = subwrd(args,9)
varUnit = subwrd(args,10)
varInst = subwrd(args,11)
varEditionRev = subwrd(args,12)

* The input variables are mprefix%' '%ssn
*                     and maveprefix%' '%ssn
*                     and oprefix%' '%ssn
*                     and oaveprefix%' '%ssn
*                     where ssn = ann djf son jja mam

titleobs_src = "$obs_title_src"
titleobs_dsc = "$obs_title_dsc"
titleobs_yrs = "$obs_title_yrs"
titleobs_att = "$obs_title_att"
titlemod     = "$mod_title"
titlemod_yrs = "$mod_title_yrs"

if ( $do_stats_file = true )
  stats_file = "$stats_file"
  rc = write(stats_file,'Var: 'varname' 'varUnit)
if( titleobs_dsc = '' )
  rc = write(stats_file,'Obs: 'titleobs_src)
else
  rc = write(stats_file,'Obs: 'titleobs_dsc)
endif
if( titleobs_att != '' )
  rc = write(stats_file,'Obs: Product Attribution: 'titleobs_att)
endif
  rc = write(stats_file,'Mod: 'titlemod)
  rc = write(stats_file,'Obs grid: 'obs_grid_type)
  rc = write(stats_file,'Mod grid: 'mod_grid_type)
  rc = write(stats_file,'Stats region: '_lon1' '_lon2' '_lat1' '_lat2)
endif

i = 1
while (i <= 5)

 'enable print PR'i

  if (i=1); ssn='ann'; endif
  if (i=2); ssn='djf'; endif
  if (i=3); ssn='jja'; endif
  if (i=4); ssn='mam'; endif
  if (i=5); ssn='son'; endif
  if (i=1); ssnUP='ANN'; endif
  if (i=2); ssnUP='DJF'; endif
  if (i=3); ssnUP='JJA'; endif
  if (i=4); ssnUP='MAM'; endif
  if (i=5); ssnUP='SON'; endif

  say '--- 'ssnUP' ---'

 'define model = 'mprefix%' '%ssn
 'define avemod = 'maveprefix%' '%ssn
 'define obs = 'oprefix%' '%ssn
 'define aveobs = 'oaveprefix%' '%ssn
 'run calStats '_lat1' '_lat2' '_lon1' '_lon2' 'model' 'obs' 'avemod' 'aveobs
  corr = subwrd(result,1)
  rmse = subwrd(result,2)
  bias = subwrd(result,3)

* top panel (model)

 'set vpage 1. 7. 7. 10.5'
  display_var(model,varname,ssnUP)
 'set annot 1 3'
 'draw title $three_back_slashes 'titlemod' 'titlemod_yrs
  plot_grid_dsc(mod)
 'run plotmeanmod avemod'
  if ( $do_stats_file = true )
    meanmod = subwrd(result,1)
    say 'meanmod = 'meanmod
  endif
 'run plotsd den1'
  if ( $do_stats_file = true )
    sdevmod = subwrd(result,1)
    say 'sdevmod = 'sdevmod
  endif

* middle panel (obs)

 'set vpage 1. 7. 3.75 7.25'
  display_var(obs,varname,ssnUP)
 'set annot 1 3'
 'draw title $three_back_slashes 'titleobs_src' 'titleobs_yrs
  plot_grid_dsc(obs)
 'run plotmeanobs aveobs'
  if ( $do_stats_file = true )
    meanobs = subwrd(result,1)
    say 'meanobs = 'meanobs
  endif
 'run plotsd den2'
  if ( $do_stats_file = true )
    sdevobs = subwrd(result,1)
    say 'sdevobs = 'sdevobs
  endif

* bottom panel (model-obs)

 'set vpage 1. 7. 0.5 4.'
 'define diff = model - obs'
  display_var(diff,varname'diff',ssnUP)
 'set annot 1 3'
  if ( strlen(''titlemod' minus 'titleobs_src'') > 40 )
    'draw title $two_back_slashes 'titlemod' minus \ 'titleobs_src
  else
    'draw title $three_back_slashes 'titlemod' minus 'titleobs_src
  endif
 'run plotStatsRegion '_lon1' '_lon2' '_lat1' '_lat2
 'run plotcorr 'corr
 'run plotbias 'bias
 'run plotrmse 'rmse
  say 'bias corr rmse = 'bias' 'corr' 'rmse

* annotation

 'set vpage 1. 7. 7.2 10.9'
 'set string 1 tc 13'
 'set strsiz 0.2 0.2'
 do_page_title = 'yes'
 if( varname = LWCRESFC )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (LWDN - LWUP) - (LWDN_clr - LWUP_clr) at SFC'
   do_page_title = 'no'
 endif
 if( varname = SWCRESFC )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (SWDN - SWUP) - (SWDN_clr - SWUP_clr) at SFC'
   do_page_title = 'no'
 endif
 if( varname = NETCRESFC )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (SWCRE + LWCRE) at SFC'
   do_page_title = 'no'
 endif
 if( varname = NETLWSFC )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (LWDN - LWUP) at SFC'
   do_page_title = 'no'
 endif
 if( varname = SWABSSFC )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (SWDN - SWUP) at SFC'
   do_page_title = 'no'
 endif
 if( varname = NETRADSFC )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (SWDN - SWUP) + (LWDN - LWUP) at SFC'
   do_page_title = 'no'
 endif
 if( varname = LWUP_DIFF )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (LWUP_TOA - LWUP_SFC)'
   do_page_title = 'no'
 endif
 if( varname = SWUP_DIFF )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (SWUP_TOA - SWUP_SFC)'
   do_page_title = 'no'
 endif
 if( varname = NETLW_DIFF )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 LWUP_TOA - (LWUP_SFC - LWDN_SFC)'
   do_page_title = 'no'
 endif
 if( varname = NETSW_DIFF )
   'draw string 4.4 5.07 'ssnUP' 'varname' 'varUnit''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (SWDN_TOA - SWUP_TOA) - (SWDN_SFC - SWUP_SFC)'
   do_page_title = 'no'
 endif
 if( varname = SWCRE_RATIO )
   'draw string 4.4 5.07 'ssnUP' 'varname''
   'set strsiz 0.15 0.15'
   'draw string 4.4 4.72 (SWCRE_SFC / SWCRE_TOA)'
   do_page_title = 'no'
 endif
 if( do_page_title = 'yes' )
   'draw string 4.4 5.0 'ssnUP' 'varname' 'varUnit''
 endif
 'set vpage 0. 8.5 0.0 0.5'
 'set string 1 c 4'
 'set strsiz 0.08 0.08'
 if( titleobs_dsc != '' & titleobs_att != '' )
   'draw string 4.25 0.41 'titleobs_dsc':'
   'draw string 4.25 0.22 Publication of these observations requires CERES attribution. See:'
   'draw string 4.25 0.07 'titleobs_att
 endif

* print plot to GrADS metafile

 'print'
 'disable print PR'i

* print plot to GrADS png file (has antialiasing
* "on" but no horizontal lines in shading)

 'gxyat -x 834 -y 1080 PR'i'.png'

* clear the frame

 'c'

  if ( $do_stats_file = true )

* write out statistics into an ASCII file

    rc = write(stats_file,''ssnUP' Mod ave = 'meanmod)
    rc = write(stats_file,''ssnUP' Mod SDev = 'sdevmod)
    rc = write(stats_file,''ssnUP' Obs ave = 'meanobs)
    rc = write(stats_file,''ssnUP' Obs SDev = 'sdevobs)
    rc = write(stats_file,''ssnUP' r(Obs, Mod) = 'corr)
    rc = write(stats_file,''ssnUP' Mod - Obs = 'bias)
    rc = write(stats_file,''ssnUP' rmse = 'rmse)

  endif

 i = i + 1
endwhile

return

function display_var(var,varname,ssnUP)
*
*  technique to shade missing (undef) values a specific color
*  from a post to the GrADS Listserv by Arindam Chakraborty on
*  Thu, 23 Jan 2003: Arindam Chakraborty <arch@CAOS.IISC.ERNET.IN>
*
 'set grads off'
 'set gxout shaded'
 'd 'var
 'q gxinfo'
  xlin = sublin(result,3)
  ylin = sublin(result,4)
  xl = subwrd(xlin,4)
  xh = subwrd(xlin,6)
  yl = subwrd(ylin,4)
  yh = subwrd(ylin,6)
 'set rgb 16 220 220 220'
 'set line 16'
 'draw recf 'xl' 'yl' 'xh' 'yh
 'run clevs_obsmod 'varname' 'ssnUP
 'd 'var
 'run cbar.gs'
return

function plot_grid_dsc(field_var)

* annotate field variable plot with a description of its grid

* extract edge of plot

  'q gxinfo'
  info=result
  xinfo=sublin(info,3)
  yinfo=sublin(info,4)
  xl=subwrd(xinfo,4) ; xr=subwrd(xinfo,6)
  yt=subwrd(yinfo,4)
  xp=xl+0.5*(xr-xl)
  yp=0.5*yt

* annotate plot

  'set string 1 c 3'
  'set strsiz 0.11 0.11'
  if( field_var = mod )
    'draw string 'xp' 'yp' '"$mod_grid_dsc"
    return
  endif
  if( field_var = obs )
    'draw string 'xp' 'yp' '"$obs_grid_dsc"
    return
  endif

return
EOF

cd $local_run_dir || exit 1

echo
echo 'GrADS script files in run work directory:'
echo $cwd
echo
ls -lt *.gs || exit 1

echo
echo 'GrADS plot_3panel_radiation.gs file:'
echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
cat plot_3panel_radiation.gs
echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
echo

foreach plot_file ( $local_plot_file_dir/*.plot )

  set variable = $plot_file:t
  set variable = $variable:r

  echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  echo '                      '$variable
  echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  echo
  echo 'working in directory:'
  echo $cwd

  set local_out_dir = $local_rad_out_rdir/$variable
  mkdir -p $local_out_dir || exit 1

  echo
  echo 'GrADS '$plot_file' file:'
  echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
  cat $plot_file
  echo '-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -'
  echo
  echo 'run GrADS with '$plot_file' file as input...'

  grads -bp < $plot_file || exit 1

  gxps -c -i PR1 -o $local_out_dir/$variable.ann.ps || exit 1
  gxps -c -i PR2 -o $local_out_dir/$variable.djf.ps || exit 1
  gxps -c -i PR3 -o $local_out_dir/$variable.jja.ps || exit 1
  gxps -c -i PR4 -o $local_out_dir/$variable.mam.ps || exit 1
  gxps -c -i PR5 -o $local_out_dir/$variable.son.ps || exit 1

  mv PR1.png $local_out_dir/$variable.ann.png || exit 1
  mv PR2.png $local_out_dir/$variable.djf.png || exit 1
  mv PR3.png $local_out_dir/$variable.jja.png || exit 1
  mv PR4.png $local_out_dir/$variable.mam.png || exit 1
  mv PR5.png $local_out_dir/$variable.son.png || exit 1

  # remove GrADS meta files

  rm -f PR? || exit 1

  # compress postscript files

  if ( ! $?NoCompress ) then
    gzip -fr $local_out_dir/$variable.???.ps || exit 1
  endif

  # move stats file to $local_out_dir/$variable.stats.txt

  if ( $do_stats_file == true ) then
    mv $stats_file $local_out_dir/$variable.stats.txt || exit 1
  endif

  # move plot file and copy plot_3panel_radiation.gs to $local_out_dir

  mv $plot_file $local_out_dir || exit 1
  cp -p plot_3panel_radiation.gs $local_out_dir || exit 1

  if ( $?local_rad_out_rdir2 ) then
    if ( $variable == 'netrad.ceres_ebaf_ed2.6' ) then
      set local_out_dir2 = $local_rad_out_rdir2/Seman.netradtoa.ceres_ebaf_ed2.6
      mkdir -p $local_out_dir2 || exit 1
      cp -p $local_out_dir/* $local_out_dir2 || exit 1
      echo
      echo 'contents of:'
      echo $local_out_dir
      echo 'copied to:'
      echo $local_out_dir2
      cd $local_out_dir2 || exit 1
      foreach file ( * )
        set to_file = `echo $file | sed s/netrad/netradtoa/g` || exit 1
        if ( $to_file != $file ) mv $file $to_file || exit 1
      end
      echo 'and "netrad" changed to "netradtoa" in names of files in:'
      echo $local_out_dir2
      cd $local_run_dir || exit 1
    else
      set local_out_dir2 = $local_rad_out_rdir2/Seman.$variable
      mkdir -p $local_out_dir2 || exit 1
      cp -p $local_out_dir/* $local_out_dir2 || exit 1
      echo
      echo 'contents of:'
      echo $local_out_dir
      echo 'copied to:'
      echo $local_out_dir2
    endif
  endif

  echo

end	# radiation variable

  if ( $?local_rad_out_rdir2 ) then
    unset local_rad_out_rdir2
  endif

end	# obs radiation type

# reference for following gcp code: Chris Golaz, 9/10/2013

echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo 'transfer output files from local to remote output directory...'
gcp -cd -r $local_out_rdir/* $out_dir || exit 1
echo '...output files transferred from local output directory:'
echo $local_out_rdir
echo 'to remote output directory:'
echo $out_dir
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

echo

echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo '...finished.'
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo

exit 0
