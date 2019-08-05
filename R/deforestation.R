#' 2000-2010 deforestation rates for the northeastern ecuadorian Amazon.
#'
#' This dataset contains the annual deforestation rate acoording FAO (FAO, 2003) from 2000-2010 for 2418 grid cells (each one of 400 ha.). Additionally, 35 variables from landscape, commodities, socioeconomic and sociocultural features are included. Sources are: Instituto Nacional de Estadisticas y Censos (INEC 2001, 2010), Sistema Nacional de Informacion (SNI, 2017), Sistema Nacional de Informacion y Gestion de Tierras Rurales e Infraestructura Tecnologica (SIGTIERRAS, 2015), National Oceanic and Atmospheric Administration (NOAA, 2019) and deforestation maps from Santos et al. 2018.
#'
#' @format A SpatialPolygonsDataFrame with 2418 rows and 37 variables:
#' \describe{
#'   \item{ID_grid}{an unique numeric identifier for each grid cell}
#'   \item{fao}{annual deforestation rate (2000-2010)}
#'   \item{A_cao}{Accessibility to coffee and cacao collection centers (hours)}
#'   \item{A_fru}{Accessibility to fruits collection centers (hours)}
#'   \item{A_mlk}{Accessibility to milk collection centers (hours)}
#'   \item{A_plm}{Accessibility to palm oil extraction facilities (hours)}
#'   \item{I_min}{Distance to mining blocks assigned between 2000 and 2010 (meters)}
#'   \item{I_ngt}{Stable nightlights trend 2000-2010 (slope)}
#'   \item{I_oil}{Distance to oil wells perforated between 2000 and 2010 (meters)}
#'   \item{B_alt}{Altitude (meters above sea level)}
#'   \item{B_rfl}{Annual rainfall (mm)}
#'   \item{B_fer}{Soil fertility (percentage organic matter)}
#'   \item{C_bsl}{Bare soil (percentage frequency)}
#'   \item{C_fra}{Fractal dimension index (unitless)}
#'   \item{C_pas}{Pasture (percentage frequency)}
#'   \item{C_sze}{Mean patch size (ha)}
#'   \item{D_adt}{Adult population (26-45 yrs)}
#'   \item{D_old}{Older adult population (45-72 yrs)}
#'   \item{D_ygr}{Young population (15-25 yrs)}
#'   \item{E_hgr}{Higher education (>13 yrs)}
#'   \item{E_ilt}{Illiterate}
#'   \item{E_pri}{Primary education (1-6 yrs)}
#'   \item{E_sec}{Secondary education (7-12 yrs)}
#'   \item{G_chf}{Chief female household}
#'   \item{G_chm}{Chief male household}
#'   \item{G_pof}{Female population}
#'   \item{G_pom}{Male population}
#'   \item{H_lar}{Large families (>5 children)}
#'   \item{H_med}{Medium families (3-5 children)}
#'   \item{H_sma}{Small families (1-2 children)}
#'   \item{L_kcw}{Speak Kichwa}
#'   \item{L_oth}{Speak other languages}
#'   \item{L_spa}{Speak Spanish}
#'   \item{L_wao}{Speak Huao Tededo}
#'   \item{W_agr}{Agricultural workers}
#'   \item{W_ind}{Industrial workers}
#'   \item{W_ser}{Service workers}
#'   ...
#' }
#' @source \url{www.sigtierras.gob.ec/; www.sni.gob.ec}
"deforestation"
