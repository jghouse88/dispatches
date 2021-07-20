import numpy as np
def f(*X):
    pmax= X[0]
    pmin= X[1]
    ramp_rate= X[2]
    min_up_time= X[3]
    min_down_time= X[4]
    marg_cst= X[5]
    no_load_cst= X[6]
    st_time_hot= X[7]
    st_time_warm= X[8]
    st_time_cold= X[9]
    st_cst_hot= X[10]
    st_cst_warm= X[11]
    st_cst_cold= X[12]
    return  0.88061901344827542281735 * pmax - 0.56366029424750108134390 * pmin - 1.4615636239468736690128 * ramp_rate + 0.67097422857155687020425E-002 * min_up_time + 0.83265121444820802687481E-001 * marg_cst - 0.96767303024194832594684E-002 * no_load_cst + 0.27139056846710812864742 * st_time_hot - 0.12161717587383390204447 * st_time_warm - 0.97593668957726600887703E-001 * st_time_cold - 0.20387078695426322227924 * st_cst_hot + 0.13577840721559553127662 * pmax**2 + 0.19820407430751235677846 * pmin**2 + 1.3291405587919613573433 * ramp_rate**2 - 0.75493721327451798058794E-002 * min_up_time**2 + 0.51645353389857852344225E-001 * marg_cst**2 + 0.52287898401867840408874E-002 * no_load_cst**2 - 0.53980859184931517802397 * st_time_hot**2 - 0.32957720613192054148755 * pmax*pmin - 1.1630748451286903044632 * pmax*ramp_rate + 0.33440330698224295102872E-001 * pmax*min_up_time + 0.43882483928852064614112E-001 * pmax*min_down_time - 0.31605173351734475173380E-001 * pmax*marg_cst - 0.10467019918504957831651 * pmax*no_load_cst - 0.32896743960266411344051 * pmax*st_time_hot + 0.17298902077138217370234 * pmax*st_time_warm - 0.70080048563294403130008E-001 * pmax*st_time_cold + 0.13291484753121424189359 * pmax*st_cst_hot + 0.87498999613996841784314 * pmin*ramp_rate - 0.18879392803491351626732E-001 * pmin*min_up_time - 0.28436935307403540112992E-001 * pmin*min_down_time - 0.61603496821709662761846E-002 * pmin*marg_cst + 0.62472037068269967163836E-001 * pmin*no_load_cst + 0.18319908316282937366104 * pmin*st_time_hot - 0.93623138258127819311127E-001 * pmin*st_time_warm + 0.40701479552985690701927E-001 * pmin*st_time_cold - 0.85772508138720199299954E-001 * pmin*st_cst_hot - 0.44754462828432227394782E-001 * ramp_rate*min_up_time - 0.56803342812735517497469E-001 * ramp_rate*min_down_time + 0.33258421281662273183422E-001 * ramp_rate*marg_cst + 0.13118692033669959728925 * ramp_rate*no_load_cst + 0.44377309522981644995809 * ramp_rate*st_time_hot - 0.22802614160439624302334 * ramp_rate*st_time_warm + 0.97866684852271887407049E-001 * ramp_rate*st_time_cold - 0.19869408263765880873208 * ramp_rate*st_cst_hot - 0.91379179676388932324071E-002 * marg_cst*no_load_cst - 0.13025162909327026722339E-001 * marg_cst*st_cst_hot + 0.12814045806303691818484E-001 * no_load_cst*st_time_hot - 0.17080339704843187226269E-001 * no_load_cst*st_time_warm - 0.68000306771793167515128E-001 * (pmax*ramp_rate)**2 + 0.50584952823051483605798E-001 * (pmax*st_time_hot)**2 + 0.88388719676192294139039E-002 * (pmax*st_time_cold)**2 - 0.53958696068499352460623E-001 * (pmin*ramp_rate)**2 - 0.79168589696822274509591E-002 * (pmin*marg_cst)**2 - 0.10937475492523286366153E-001 * (pmin*st_time_hot)**2 + 0.94742645297225414752207E-002 * (ramp_rate*min_up_time)**2 - 0.12427742545633925497217E-001 * (ramp_rate*marg_cst)**2 - 0.95716927979186353786512E-002 * (ramp_rate*no_load_cst)**2 - 0.20631356514346319008801E-001 * (ramp_rate*st_time_hot)**2 + 0.34835146280662337980871E-001 * (ramp_rate*st_time_warm)**2 + 0.24431765546675789091413E-001 * (ramp_rate*st_time_cold)**2 - 0.59964112558759212479043E-001 * (ramp_rate*st_cst_hot)**2 - 0.72693114375904895399505E-002 * (min_up_time*marg_cst)**2 + 0.10564299695607940257625E-001 * (marg_cst*st_cst_hot)**2 
