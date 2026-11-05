
using WannierInterActionBase

xmldata = read("./pwscf.save/data-file-schema.xml", QE_xml_new; period = [1, 1, 1])
kline = Kline(xmldata["kpoints"], xmldata["cell"].lattice)

write("./band.dat", xmldata["eigenvalues"], BAND_dat; kline)
# write("./band_up.dat", xmldata["eigenvalues_up"], BAND_dat; kline)
# write("./band_dn.dat", xmldata["eigenvalues_dn"], BAND_dat; kline)
