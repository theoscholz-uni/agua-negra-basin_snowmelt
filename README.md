# agua-negra-basin_snowmelt
R workflow for snowmelt modelling in the Agua Negra catchment, integrating watershed delineation, ERA5 climate data, camera-derived snow cover, potential direct solar radiation, Monte Carlo calibration, and comparison with VIIRS snow-cover products. NOTE: Notations in the script are written in german!

#Bachelor thesis - University of Bonn

R code for the Agua Negra Basin workflow of my bachelor's thesis

References:

ANDERSON, E. (2006): Snow Accumulation and Ablation Model – SNOW-17. [Online-Publikation]. [Zugriff am: 12.04.2026]. Verfügbar unter: https://www.weather.gov/media/owp/oh/hrl/docs/22snow17.pdf

BEVEN, K. & BINLEY, A. (1992): The future of distributed models: Model calibration and uncertainty prediction. In: Hydrological Processes 6, 3. S. 279-298.

DALLA TORRE, D., DI MARCO, N., MENAPACE, A., AVESANI, D., RIGHETTI, M. & MAJONE, B. (2024): Suitability of ERA5-Land reanalysis dataset for hydrological modelling in the Alpine region. In: Journal of Hydrology: Regional Studies 52. DOI: 10.1016/j.ejrh.2024.101718.

ESSERY, R. & POMEROY, J. (2004): Implications of spatial distributions of snow mass and melt rate for snow-cover depletion: theoretical considerations. In: Annals of Glaciology 38. S. 261–265.

FEDOROV, R., CAMERADA, A., FRATERNALI, P. & TAGLIASACCHI, M. (2015): Estimating snow cover from publicly available images. [Online-Publikation]. [Zugriff am: 06.04.2026]. Verfügbar unter: https://arxiv.org/pdf/1508.01055

HARPOLD, A. A., KAPLAN, M. L., KLOS, P. Z., LINK, T., MCNAMARA, J. P., RAJAGOPAL, S., SCHUMER, R. & STEELE, C. M. (2017): Rain or snow: hydrologic processes, observations, prediction, and research needs. In: Hydrology and Earth System Sciences 21. S. 1-22.

HAUGHTON, N., ABRAMOWITZ, G., PITMAN, A. J., OR, D., BEST, M. J., JOHNSON, H. R., BALSAMO, G., BOONE, A., CUNTZ, M., DECHARME, B., DIRMEYER, P. A., DONG, J., EK, M., GUO, Z., HAVERD, V., VAN DEN HURK, B. J. J., NEARING, G. S., PAK, B., SANTANELLO JR., J. A., STEVENS, L. E. & VUICHARD, N. (2016): The Plumbing of Land Surface Models: Is Poor Performance a Result of Methodology or Data Quality? In: Journal of Hydrometeorology 17, 6. S. 1705–1720.

HOCK, R. (1999): A distributed temperature-index ice- and snowmelt model including potential direct solar radiation. In: Journal of Glaciology 45, 149. S. 101–111.

HÄRER, S., BERNHARDT, M. & SCHULZ, K. (2016): PRACTISE – Photo Rectification And ClassificaTIon SoftwarE (V.2.1). In: Geoscientific Model Development 9, 1. S. 307-321.

IQBAL, M. (1983): An Introduction to Solar Radiation. (Academic Press) Toronto.

ISMAIL, M. F., BOGACKI, W., DISSE, M., SCHÄFER, M. & KIRSCHBAUER, L. (2023): Estimating degreeday factors of snow based on energy flux components. In: The Cryosphere 17. S. 211–231.

JENÍČEK, M., BEITLEROVÁ, H., HASA, M., KUČEROVÁ, D., PEVNÁ, H. & PODZIMEK, S. (2012): Modelling snow accumulation and snowmelt runoff – present approaches and results. In: AUC Geographica 47, 2. S. 15–24.

JENSON, S. K. & DOMINGUE, J. O. (1988): Extracting Topographic Structure from Digital Elevation Data for Geographic Information System Analysis. In: Photogrammetric Engineering and Remote Sensing 54, 11. S. 1593-1600.

KASTEN, F. & YOUNG, A. T. (1989): Revised optical air mass tables and approximation formula. In: Applied Optics 28, 22. S. 4735-4738.

KNOBEN, W. J. M., FREER, J. E. & WOODS, R. A. (2019): Technical note: Inherent benchmark or not? Comparing Nash–Sutcliffe and Kling–Gupta efficiency scores. In: Hydrology and Earth System Sciences 23. S. 4323-4331.

LECOMTE, K. L., MILANA, J. P., FORMICA, S. M. & DEPETRIS, P. J. (2008): Hydrochemical appraisal of ice- and rock-glacier meltwater in the hyperarid Agua Negra drainage basin, Andes of Argentina. In: Hydrological Processes 22, 13. S. 2180-2195.

MOTEVALLI, M., ZADBAR, A., ELYASI, E. & JALAAL, M. (2015): Using Monte-Carlo approach for analysis of quantitative and qualitative operation of reservoirs system with regard to the inflow uncertainty. In: Journal of African Earth Sciences 105. S. 1–16.

MUÑOZ-SABATER, J., DUTRA, E., AGUSTÍ-PANAREDA, A., ALBERGEL, C., ARDUINI, G., BALSAMO, G., BOUSSETTA, S., CHOULGA, M., HARRIGAN, S., HERSBACH, H., MARTENS, B., MIRALLES, D. G., PILES, M., RODRÍGUEZ-FERNÁNDEZ, N. J., ZSOTER, E., BUONTEMPO, C. & THÉPAUT, J.-N. (2021): ERA5- Land: a state-of-the-art global reanalysis dataset for land applications. In: Earth System Science Data 13. S. 4349-4383.

O’CALLAGHAN, J. F. & MARK, D. M. (1984): The Extraction of Drainage Networks from Digital Elevation Data. In: Computer Vision, Graphics, and Image Processing 28. S. 323-344.

OUGAHI, J. H. & ROWAN, J. S. (2024): Combining Hydrological Models and Remote Sensing to Characterize Snowpack Dynamics in High Mountains. In: Remote Sensing 16, 264. DOI: 10.3390/rs16020264.

PITTE, P., MASIOKAS, M., GARGANTINI, H., RUIZ, L., BERTHIER, E., FERRI HIDALGO, L., ZALAZAR, L., DUSSAILLANT, I., VIALE, M., ZORZUT, V., CORVALÁN, E., SCARPA, J. P., COSTA, G. & VILLALBA, R. (2022): Recent mass-balance changes of Agua Negra glacier (30°S) in the Desert Andes of Argentina. In: Journal of Glaciology 68, 272. S. 1197-1209.

PORTENIER, C., HÜSLER, F., HÄRER, S. & WUNDERLE, S. (2020): Towards a webcam-based snow cover monitoring network: methodology and evaluation. In: The Cryosphere 14. S. 1409–1423.

RIGGS, G. A., HALL, D. K. & ROMÁN, M. O. (2017): Overview of NASA’s MODIS and Visible Infrared Imaging Radiometer Suite (VIIRS) snow-cover Earth System Data Records. In: Earth System Science Data 9. S. 765-777.

SALVATORI, R., PLINI, P., GIUSTO, M., VALT, M., SALZANO, R., MONTAGNOLI, M., CAGNATI, A., CREPAZ, G. & SIGISMONDI, D. (2011): Snow cover monitoring with images from digital camera systems. In: Italian Journal of Remote Sensing 43, 2. S. 137–145.

SCHROTT, L. (1991): Global solar radiation, soil temperature and permafrost in the Central Andes, Argentina: A progress report. In: Permafrost and Periglacial Processes 2. S. 59-66.

SETIYOWATI, Y. A., HARISUSENO, D. & SAJALI, M. A. (2025): Comparison of correlation, PBIAS and RSR between monthly, daily, and hourly GPM rainfall data. In: Jurnal Penelitian Pendidikan IPA 11, 6. S. 615–624.

SHIELDS, M. D. & ZHANG, J. (2016): The generalization of Latin hypercube sampling. In: Reliability Engineering and System Safety 148. S. 96–108.

STRANG, G. (2023): Calculus. [Online-Publikation]. [Zugriff am: 08.04.2026]. Verfügbar unter: https://ocw.mit.edu/courses/res-18-001-calculus-fall-2023/

SZOKOLAY, S. V. (2007): Solar Geometry. [Online-Publikation]. [Zugriff am: 06.04.2023]. Verfügbar unter: https://www.plea-arch.org/wp-content/uploads/PLEA-NOTE-1-SOLAR-GEOMETRY.pdf

ZHANG, F., ZHANG, H., HAGEN, S. C., YE, M., WANG, D., GUI, D., ZENG, C., TIAN, L. & LIU, J. (2015): Snow cover and runoff modelling in a high mountain catchment with scarce data: effects of temperature and precipitation parameters. In: Hydrological Processes 29. S. 52-65.

ZHANG, T., STACKHOUSE, P. W., MACPHERSON, B. & MIKOVITZ, J. C. (2021): A solar azimuth formula that renders circumstantial treatment unnecessary without compromising mathematical rigor: Mathematical setup, application and extension of a formula based on the subsolar point and atan2 function. In: Renewable Energy 172. S. 1333-1340.
