# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 18:29:19 2015

@author: k
"""

def cleanms1(filename):
    """
    given a txt of orbitrap mass spect file.  for each scan, if it is MS scan, not MS2 scan, then delete the DataPeaks part
    """
    fo = open(filename,"r")
    line = fo.readline()
    fo2 = open(filename+"new","w")
    while "ScanHeader # " not in line:
        if line != "\n":
            fo2.write(line)
        line = fo.readline()
    scan = line
    line = fo.readline()
    def cleanscan(scan):
        """
        scan is all information from one scan
        if it is MS scan, delete the DataPeaks part. else output the whole part
        """
        if "Scan Type, MS Scan" in scan:
            return scan.split("DataPeaks\n")[0]
        else:
            return scan
    allresult =[]
    while line:
        if "ScanHeader # " not in line:
            if line != "\n":
                scan += line
            line = fo.readline()
        else:
            scanout = cleanscan(scan)
            allresult.append(scanout)
            scan = line
            line = fo.readline()
    for ele in allresult:
        fo2.write(ele)
    fo2.close()
    fo.close()

folder = "E:\\store_for_D\\ManducaPeptidesProteome\\RawTxt\\"
file = "1380_Control_mSIM2.txtnew"
cleanms1(folder+file)


file = "1380_Control_mSIM3.txtnew"
cleanms1(folder+file)
file = "1380_Control_mTP1.txtnew"
cleanms1(folder+file)
file = "1380_Control_mTP3.txtnew"
cleanms1(folder+file)
file = "1380_Induced_m100.txtnew"
cleanms1(folder+file)
file = "1380_Induced_mTP1.txtnew"
cleanms1(folder+file)
file = "1380_Induced_mTP3.txtnew"
cleanms1(folder+file)

file = "1380_Induced_mSIM2.txtnew"
cleanms1(folder+file)
file = "1380_Induced_mSIM3.txtnew"
cleanms1(folder+file)
file2 ="1380_Control_m100min.txtnew"
cleanms1(folder+file2)



#20160826
#synthetic peptide, data processing.
filename = r"E:\Lab\works\20160516MSMS_sytheticPeptides\1437_XC_peptides_1.ms1"
outfolder = "E:\\Lab\\works\\20160516MSMS_sytheticPeptides\\scans\\"
import ms1Class
ms1Class.ms1ToIndividualScans(filename, outfolder)

#20160828 filter the ms1 file
def filterms1Files():
    import ms1Class
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.ms1"
    ms1Class.ms1CleanUp(filename)
    
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.ms1"
    ms1Class.ms1CleanUp(filename)
    
def readms1fileToDC():
    import ms1Class
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.ms1.clear"
    dcmzs, dcmzi = ms1Class.readms1ToMzsMzi(filename)
    
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.ms1.clear"
    kbmzs, kbmzi = ms1Class.readms1ToMzsMzi(filename)


def massMs1Dealwith20160829():
    """
    deal with the data we got from mass spec
    """
    peptidesfile = r"E:\Lab\works\20160516MSMS_sytheticPeptides\analysis\10pep.fasta"
    from Bio import SeqIO
    peptides = list(SeqIO.parse(open(peptidesfile),'fasta'))
    import chemicalAllPossibleMsms
    dc_pep = {}
    for ele in peptides:
        dc_pep[ele.id] = set(chemicalAllPossibleMsms.modiPeptide(str(ele.seq)) + \
            chemicalAllPossibleMsms.modiPeptide(str(ele.seq), dcaaModi = chemicalAllPossibleMsms.DCaaModi2))
    
    dc_mz ={} #keep three digit after decimal point
    charge_max = 5
    for ele in dc_pep:
        for elepep in dc_pep[ele]:
            dc_mz[elepep] = [round(chemicalAllPossibleMsms.getPepmassWithModi(elepep, charge=num),3) for num in range(1, charge_max+1)]
    
    dc_mz_neutron = {} #key (peptide, charge), value:[mz *4] with 0 to 3 neutrons
    neutron_max = 4
    neutron_w = 1.00727
    for ele in dc_pep:
        for elepep in dc_pep[ele]:
            for charge in range(1, charge_max+1):
                mz = chemicalAllPossibleMsms.getPepmassWithModi(elepep, charge=charge)
                mzs = [round(mz + neutron_num * neutron_w / charge,3) for neutron_num in range(neutron_max +1)]
                dc_mz_neutron[(elepep, charge)] = mzs
    
    mz_interest = set()
    for ele in dc_mz_neutron.values():
        mz_interest.update(ele)
    mz_max = 1400
    mz_min = 300
    mz_interest = [ele for ele in mz_interest if ele <=mz_max and ele >= mz_min]
    mz_interest.sort()
    
    #get numpy array for synthetic peptide ms1 result
    import ms1Class
    mzi_min = 10000
    mzi_min_len = 3 # at least three adjacent mzis not zero
    mzis_min = 100000
    ms1file1437 = r"E:\Lab\works\20160516MSMS_sytheticPeptides\1437_XC_peptides_1.ms1"
    npmzs1437, npmzi1437, scan_info1437 = ms1Class.readms1ToMzsMzi(ms1file1437, mzi_min, True)
    
    dc_mzscan = {}
    if mzi_min_len<=1:
        mzi_min_len = 1
    for ele in mz_interest:
        templsmzi = ms1Class.targetmzIntensityFromNpmzsmzi(ele, npmzs1437, npmzi1437)
        #clean dc_mzscan if max mzi is less than 100,000
        if max(templsmzi) > mzis_min:
            for num in range(len(templsmzi)-mzi_min_len+1):
                over0 = True
                for num2 in range(num, num + mzi_min_len):
                    if templsmzi[num2] == 0:
                        over0 = False
                        break
                if over0:
                    dc_mzscan[ele] = templsmzi
                    break
    
    #find good peptides
    #find good groups in dc_mz_neutron:
    dc_mz_neutron_good = {}
    for key in dc_mz_neutron:
        good = True
        for mz in dc_mz_neutron[key]:
            if mz not in dc_mzscan:
                good = False
                break
        if good:
            dc_mz_neutron_good[key] = dc_mz_neutron[key]
    
    peptides_good = [ele[0] for ele in list(dc_mz_neutron_good.keys())]
    peptides_good = set(peptides_good)
    dc_pep_final = {}
    mz_track = []
    for name in dc_pep:
        dc_pep_final[name] = []
        for peptide in dc_pep[name]:
            if peptide in peptides_good:
                for charge in range(1, charge_max+1):
                    if (peptide, charge) in dc_mz_neutron_good:
                        mz = dc_mz_neutron_good[(peptide, charge)][0]
                        if mz not in mz_track:
                            mz_track.append(mz)
                            dc_pep_final[name].append((peptide, charge))
    
    #plot for each peptide
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('20160829syntheticPeptidesMassSpec2.pdf')
    plt.rcParams["figure.figsize"] = [20.0,10]
    cm = plt.get_cmap('nipy_spectral')
    mzused = []
    for peptide in peptides:
        name = peptide.id
        x = [ele[1] for ele in scan_info1437]
        NUM_COLORS = len(dc_pep_final[name])
        for num in range(len(dc_pep_final[name])):
            mz = dc_mz_neutron_good[dc_pep_final[name][num]][0]
            if mz in mzused:
                continue
            mzused.append(mz)
            y = dc_mzscan[mz]
            total_intensity = sum(y)
            label = '%30s'%dc_pep_final[name][num][0] +' %5d'%dc_pep_final[name][num][1] + ' %9.3f'%mz + ' %10.2e'%total_intensity
#            print(label)
            colorcode = num/NUM_COLORS
            colorcode = colorcode - int(colorcode)
            plt.plot(x,y,label = label,color = cm(colorcode))
            plt.legend(loc = 0)
        plt.xlim(55,90)
        plt.title(name + ' seq charge m/z total intensity')
        plt.tight_layout()
        plt.xlabel('retention time')
        plt.ylabel('m/z intensity')
        pp.savefig()
        plt.close()
    pp.close()
    
    def mspeptidefinder(peptidesfile, ms1file, outpdf_file, charge_max = 5, neutron_max = 3):
        """
        change the code above to make it working for different files.
        """            
        from Bio import SeqIO
        peptides = list(SeqIO.parse(open(peptidesfile),'fasta'))
        import chemicalAllPossibleMsms
        dc_pep = {}
        for ele in peptides:
            dc_pep[ele.id] = set(chemicalAllPossibleMsms.modiPeptide(str(ele.seq)) + \
                chemicalAllPossibleMsms.modiPeptide(str(ele.seq), dcaaModi = chemicalAllPossibleMsms.DCaaModi2))
        
        dc_mz ={} #keep three digit after decimal point
        charge_max = charge_max
        for ele in dc_pep:
            for elepep in dc_pep[ele]:
                dc_mz[elepep] = [round(chemicalAllPossibleMsms.getPepmassWithModi(elepep, charge=num),3) for num in range(1, charge_max+1)]
        
        dc_mz_neutron = {} #key (peptide, charge), value:[mz *4] with 0 to 3 neutrons
        neutron_max = neutron_max
        neutron_w = 1.00727
        for ele in dc_pep:
            for elepep in dc_pep[ele]:
                for charge in range(1, charge_max+1):
                    mz = chemicalAllPossibleMsms.getPepmassWithModi(elepep, charge=charge)
                    mzs = [round(mz + neutron_num * neutron_w / charge,3) for neutron_num in range(neutron_max +1)]
                    dc_mz_neutron[(elepep, charge)] = mzs
        
        mz_interest = set()
        for ele in dc_mz_neutron.values():
            mz_interest.update(ele)
        mz_max = 1500
        mz_min = 300
        mz_interest = [ele for ele in mz_interest if ele <=mz_max and ele >= mz_min]
        mz_interest.sort()
        
        #get numpy array for synthetic peptide ms1 result
        import ms1Class
        mzi_min = 10000
        mzi_min_len = 3 # at least three adjacent mzis not zero
        mzis_min = 100000
        npmzs1437, npmzi1437, scan_info1437 = ms1Class.readms1ToMzsMzi(ms1file, mzi_min, True)
        
        dc_mzscan = {}
        if mzi_min_len<=1:
            mzi_min_len = 1
        for ele in mz_interest:
            templsmzi = ms1Class.targetmzIntensityFromNpmzsmzi(ele, npmzs1437, npmzi1437, 0.003, 3)
            #clean dc_mzscan if max mzi is less than 100,000
            if max(templsmzi) > mzis_min:
                for num in range(len(templsmzi)-mzi_min_len+1):
                    over0 = True
                    for num2 in range(num, num + mzi_min_len):
                        if templsmzi[num2] == 0:
                            over0 = False
                            break
                    if over0:
                        dc_mzscan[ele] = templsmzi
                        break
        
        #find good peptides
        #find good groups in dc_mz_neutron:
        dc_mz_neutron_good = {}
        for key in dc_mz_neutron:
            good = True
            for mz in dc_mz_neutron[key]:
                if mz not in dc_mzscan:
                    good = False
                    break
            if good:
                dc_mz_neutron_good[key] = dc_mz_neutron[key]
        
        peptides_good = [ele[0] for ele in list(dc_mz_neutron_good.keys())]
        peptides_good = set(peptides_good)
        dc_pep_final = {}
        mz_track = []
        for name in dc_pep:
            dc_pep_final[name] = []
            for peptide in dc_pep[name]:
                if peptide in peptides_good:
                    for charge in range(1, charge_max+1):
                        if (peptide, charge) in dc_mz_neutron_good:
                            mz = dc_mz_neutron_good[(peptide, charge)][0]
                            if mz not in mz_track:
                                mz_track.append(mz)
                                dc_pep_final[name].append((peptide, charge))
        
        #plot for each peptide
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(outpdf_file)
        plt.rcParams["figure.figsize"] = [20.0,10]
        cm = plt.get_cmap('nipy_spectral')
        mzused = []
        for peptide in peptides:
            name = peptide.id
            x = [ele[1] for ele in scan_info1437]
            NUM_COLORS = len(dc_pep_final[name])
            for num in range(len(dc_pep_final[name])):
                mz = dc_mz_neutron_good[dc_pep_final[name][num]][0]
                if mz in mzused:
                    continue
                mzused.append(mz)
                y = dc_mzscan[mz]
                total_intensity = sum(y)
                label = '%30s'%dc_pep_final[name][num][0] +' %5d'%dc_pep_final[name][num][1] + ' %9.3f'%mz + ' %10.2e'%total_intensity
                print(label)
                colorcode = num/NUM_COLORS
                colorcode = colorcode - int(colorcode)
                plt.plot(x,y,label = label,color = cm(colorcode))
                plt.legend(loc = 0)
            plt.xlim(20,120)
            plt.title(name + ' seq charge m/z total intensity')
    #        plt.tight_layout()
            plt.xlabel('retention time')
            plt.ylabel('m/z intensity')
            pp.savefig()
            plt.close()
        pp.close()
        
        #plot according to peptide sequence
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(outpdf_file)
        plt.rcParams["figure.figsize"] = [20.0,10]
        cm = plt.get_cmap('nipy_spectral')
        mzused = []
        dc_pep_pepBased = {}
        for ele in dc_pep_final:
            dc_pep_pepBased[ele] = {}
            for ele2 in dc_pep_final[ele]:
                if ele2[0] not in dc_pep_pepBased[ele]:
                    dc_pep_pepBased[ele][ele2[0]] = []
                dc_pep_pepBased[ele][ele2[0]].append(ele2)
        for peptide in peptides:
            name = peptide.id
            x = [ele[1] for ele in scan_info1437]
            for ele2 in dc_pep_pepBased[name]:
                keys = dc_pep_pepBased[name][ele2]    
                NUM_COLORS = len(keys)
                for num in range(len(keys)):
                    mz = dc_mz_neutron_good[keys[num]][0]
                    if mz in mzused:
                        continue
                    mzused.append(mz)
                    y = dc_mzscan[mz]
                    total_intensity = sum(y)
                    label = '%30s'%keys[num][0] +' %5d'%keys[num][1] + ' %9.3f'%mz + ' %10.2e'%total_intensity
#                    print(label)
                    colorcode = num/NUM_COLORS
                    colorcode = colorcode - int(colorcode)
                    plt.plot(x,y,label = label,color = cm(colorcode))
                    plt.legend(loc = 0)
                plt.xlim(20,120)
                plt.title(name + ' seq charge m/z total intensity')
        #        plt.tight_layout()
                plt.xlabel('retention time')
                plt.ylabel('m/z intensity')
                pp.savefig()
                plt.close()
        pp.close()
    
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\7pep.txt"
    ms1file = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.ms1.clear"
    outpdf_file = r'20160830MsBarLessthan10kdMS_3.pdf'
    charge_max = 5
    neutron_max = 4
    mspeptidefinder(peptidesfile,ms1file,outpdf_file)
    
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\seq4.fasta"
    ms1file = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.ms1.clear"
    outpdf_file = r'20160830AgKrishaHemo2.pdf'
    charge_max = 5
    neutron_max = 4
    mspeptidefinder(peptidesfile,ms1file,outpdf_file)

    peptidesfile = r'E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\Yan\target.txt'
    ms1file = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\Yan\1352_HP1.ms1"
    outpdf_file = r'20160830Yan.pdf'
    mspeptidefinder(peptidesfile,ms1file,outpdf_file,5,0)
    
    
def ms1dealWithNewMethod20160901():
    peptidesfile = r"E:\Lab\works\20160516MSMS_sytheticPeptides\analysis\10pep.fasta"
    from Bio import SeqIO
    peptides = list(SeqIO.parse(open(peptidesfile),'fasta'))
    
    import chemicalAllPossibleMsms
    dc_pep = {}
    for ele in peptides:
        dc_pep[ele.id] = set(chemicalAllPossibleMsms.modiPeptide(str(ele.seq), dcaaModi = {'C':'d'},random = False))
    
    import ms1Class
    mzi_min = 10000
    mzi_min_len = 3 # at least three adjacent mzis not zero
    ms1file = r"E:\Lab\works\20160516MSMS_sytheticPeptides\1437_XC_peptides_1.ms1mgflike"
    npmzs, npmzi, scan_info = ms1Class.readms1ToMzsMzi(ms1file, mzi_min, True)
    
    mz_min = 300
    mz_max = 1500
    charge_max =5
    neutron_max = 3
    dc_mzis = {}
    dc_mz = {}
    for name in dc_pep:
        dc_mzis[name] = {}
        dc_mz[name] = {}
        for peptide in dc_pep[name]:
            dc_mz[name][peptide] = {}
            for charge in range(1, charge_max+1):
                mz = chemicalAllPossibleMsms.getPepmassWithModi(peptide,charge)
                dc_mz[name][peptide][charge] = mz
                if mz >= mz_min and mz <= mz_max:
                    dc_mzis[name][charge] = ms1Class.targetPepmzIntensityFromNpmzsmzi(peptide,npmzs, npmzi,charge, neutron_max,0.003, mzi_min_len)
    
    outpdf_file = '20160901syntheticpeptide_HPLC.pdf'
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(outpdf_file)
    plt.rcParams["figure.figsize"] = [20.0,10]
    cm = plt.get_cmap('nipy_spectral')
    x = [ele[1] for ele in scan_info]
    for ele in peptides:
        name = ele.id
        peptide = dc_pep[name]
        NUM_COLORS = len(dc_mzis[name])
        num = 0
        for charge in dc_mzis[name]:
            y = dc_mzis[name][charge]
            total_intensity = sum(y)
            mz = chemicalAllPossibleMsms.getPepmassWithModi(peptide,charge)
            label = '%30s'%list(dc_pep[name])[0] +' %5d'%charge + ' %9.3f'%mz + ' %10.2e'%total_intensity
            if total_intensity > 0:
                colorcode = num/NUM_COLORS
                plt.plot(x,y,label = label,color = cm(colorcode))
                plt.legend(loc = 0)
                num +=1
        plt.xlim(55,90)
        plt.title(name + ' seq charge m/z total intensity')
#        plt.tight_layout()
        plt.xlabel('retention time')
        plt.ylabel('m/z intensity')
        pp.savefig()
        plt.close()
    pp.close()
    
    def mspeptidefinder(peptidesfile, ms1file, outpdf_file, dcaaModi,charge_max = 5,neutron_max = 3):
#        peptidesfile = r"E:\Lab\works\20160516MSMS_sytheticPeptides\analysis\10pep.fasta"
#        dcaaModi = {'C':'d'}
        from Bio import SeqIO
        peptides = list(SeqIO.parse(open(peptidesfile),'fasta'))
        
        import chemicalAllPossibleMsms
        dc_pep = {}
        for ele in peptides:
            dc_pep[ele.id] = set(chemicalAllPossibleMsms.modiPeptide(str(ele.seq), dcaaModi = dcaaModi,random = False))
        
        import ms1Class
        mzi_min = 10000
        mzi_min_len = 3 # at least three adjacent mzis not zero
#        ms1file = r"E:\Lab\works\20160516MSMS_sytheticPeptides\1437_XC_peptides_1.ms1mgflike"
        npmzs, npmzi, scan_info = ms1Class.readms1ToMzsMzi(ms1file, mzi_min, True)
        
        mz_min = 300
        mz_max = 1500
#        charge_max =5
#        neutron_max = 3
        dc_mzis = {}
        for name in dc_pep:
            dc_mzis[name] = {}
            for peptide in dc_pep[name]:
                for charge in range(1, charge_max+1):
                    mz = chemicalAllPossibleMsms.getPepmassWithModi(peptide,charge)
                    if mz >= mz_min and mz <= mz_max:
                        dc_mzis[name][charge] = ms1Class.targetPepmzIntensityFromNpmzsmzi(peptide,npmzs, npmzi,charge, neutron_max,0.003, mzi_min_len)
        
#        outpdf_file = '20160901syntheticpeptide_HPLC.pdf'
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(outpdf_file)
        plt.rcParams["figure.figsize"] = [20.0,10]
        cm = plt.get_cmap('nipy_spectral')
        x = [ele[1] for ele in scan_info]
        for ele in peptides:
            name = ele.id
            peptide = dc_pep[name]
            NUM_COLORS = len(dc_mzis[name])
            num = 0
            for charge in dc_mzis[name]:
                y = dc_mzis[name][charge]
                total_intensity = sum(y)
                mz = chemicalAllPossibleMsms.getPepmassWithModi(peptide,charge)
                label = '%30s'%list(dc_pep[name])[0] +' %5d'%charge + ' %9.3f'%mz + ' %10.2e'%total_intensity
                if total_intensity > 0:
                    colorcode = num/NUM_COLORS
                    plt.plot(x,y,label = label,color = cm(colorcode))
                    plt.legend(loc = 0)
                    num +=1
            plt.xlim(20,120)
            plt.title(name + ' seq charge m/z total intensity')
    #        plt.tight_layout()
            plt.xlabel('retention time')
            plt.ylabel('m/z intensity')
            pp.savefig()
            plt.close()
        pp.close()
        
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\7pep.txt"
    ms1file = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.ms1.clear"
    outpdf_file = r'20160901MsBarLessthan10kdMS_3.pdf'
    charge_max = 5
    neutron_max = 3
    dcaaModi = {'C':'d'}
    mspeptidefinder(peptidesfile,ms1file,outpdf_file,dcaaModi, charge_max,neutron_max)
    
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\seq4.fasta"
    ms1file = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.ms1.clear"
    outpdf_file = r'20160901AgKrishaHemo2.pdf'
    charge_max = 5
    neutron_max = 3
    dcaaModi = {'C':'c'}
    mspeptidefinder(peptidesfile,ms1file,outpdf_file,dcaaModi, charge_max,neutron_max)


def cleanMS1toPeaks20160907():
    import ms1Class
    filename = r"E:\Lab\works\20160516MSMS_sytheticPeptides\1437_XC_peptides_1.ms1"
    ms1Class.ms1ToPeaks(filename,None,10000)
    
    import ms1Class
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.ms1"
    ms1Class.ms1ToPeaks(filename,None,100000)
    
    import ms1Class
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.ms1.clear"
    ms1Class.ms1ToPeaks(filename,None,100000)

    import ms1Class
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_highrange.ms1"
    ms1Class.ms1ToPeaks(filename,None,100000)

    import ms1Class
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_highrange.ms1"
    ms1Class.ms1ToPeaks(filename,filename+'2Peak',100000)
    
    import ms1Class
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_highrange.ms1"
    ms1Class.ms1ToPeaks(filename,filename+'2Peak',100000)


def ms1dealWithNeutron20160908():
    
    
    def ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi = {'C':'d'}, random = False):
        #read in the fasta file    
#        peptidesfile = r"E:\Lab\works\20160516MSMS_sytheticPeptides\analysis\10pep.fasta"
        from Bio import SeqIO
        peptides = list(SeqIO.parse(open(peptidesfile),'fasta'))
        #get all possible target sequence with modification for each peptide name
        import chemicalAllPossibleMsms
        dc_pep = {}
        for ele in peptides:
            dc_pep[ele.id] = set(chemicalAllPossibleMsms.modiPeptide(str(ele.seq), dcaaModi = dcaaModi,random = random))
        
        #clean dc_pep for random modification
        for ele in dc_pep:
            dc_pep[ele] = [i for i in dc_pep[ele] if i.count('c') == 2 or i.count('c') == 0]
        
        neutron_max = 4
        charge_max = 5
        dc_mz_neutron = {}
        for name in dc_pep:
            dc_mz_neutron[name] = {}
            for peptide in dc_pep[name]:
                dc_mz_neutron[name][peptide] = {}
                for charge in range(1, 1 + charge_max +1):
                    dc_mz_neutron[name][peptide][charge] = [chemicalAllPossibleMsms.getPepmassWithModi(peptide, charge, 'M', n) for n in range(neutron_max+1)]
                        
        #readin the ms1Peak file
        import ms1Class
#        ms1Peak = r"E:\Lab\works\20160516MSMS_sytheticPeptides\1437_XC_peptides_1.ms1.Peak"
        npmzs, npmzi, scan_info = ms1Class.readms1Peaks2mzsmzi(ms1Peak)
        
        dc_mzscan = {}
        mz_min = 300
        mz_max = 2000
        peak_min = 100000
        for name in dc_mz_neutron:
            for peptide in dc_mz_neutron[name]:
                for charge in dc_mz_neutron[name][peptide]:
                    for neutron in range(neutron_max):
                        mz = dc_mz_neutron[name][peptide][charge][neutron][0]
                        if mz > mz_min and mz < mz_max:
                            lsmzi = ms1Class.targetmzIntensityFromNpmzsmzi(mz,npmzs, npmzi,0.01)
                            if max(lsmzi) > peak_min:
                                dc_mzscan[(name, peptide,charge,neutron)] = lsmzi
    
            
        #plot for each peptide
#        outfile = '20160908syntheticPeptidesMassSpec2.pdf'
        import matplotlib.pyplot as plt
#        import numpy as np
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(outfile)
        plt.rcParams["figure.figsize"] = [20.0,10]
        cm = plt.get_cmap('gist_rainbow')
        x = [ele[1] for ele in scan_info]
        for name in peptides:
            name = name.id
            for peptide in dc_mz_neutron[name]:
                for charge in dc_mz_neutron[name][peptide]:
                    ploted = False
                    for neutron in range(neutron_max):
                        mz = dc_mz_neutron[name][peptide][charge][neutron][0]
                        percent = dc_mz_neutron[name][peptide][charge][neutron][1]
                        key = (name, peptide,charge,neutron)
                        if key in dc_mzscan:
                            ploted = True
                            y = dc_mzscan[key]
                            label = ' charge:' +str(charge) + \
                            ' neutron' + str(neutron) +' mz:%.3f'%mz + ' ratio:%.3f' %percent
                            plt.plot(x,y,label = label, color = cm(neutron / 4),linewidth = 0.3)
                    if ploted:
                        plt.title(name +' ' + peptide)
                        plt.xlim(20,120)
                        plt.tight_layout()
                        plt.legend(loc = 0)
                        pp.savefig()
                        plt.close()
        pp.close()
        
        
        #for Bar10KD
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(outfile)
        plt.rcParams["figure.figsize"] = [20.0,10]
        cm = plt.get_cmap('gist_rainbow')
        x = [ele[1] for ele in scan_info]
        for name in peptides:
            name = name.id
            for peptide in dc_mz_neutron[name]:
                for charge in dc_mz_neutron[name][peptide]:
                    ploted = False
                    for neutron in range(neutron_max):
                        mz = dc_mz_neutron[name][peptide][charge][neutron][0]
                        percent = dc_mz_neutron[name][peptide][charge][neutron][1]
                        key = (name, peptide,charge,neutron)
                        if key in dc_mzscan:
                            y = dc_mzscan[key]
                            label = ' charge:' +str(charge) + \
                            ' neutron' + str(neutron) +' mz:%.3f'%mz + ' ratio:%.3f' %percent
                            if ((name == 'SRP1' or name == 'SRP4' or name == 'uENF1' and peptide.count('o') == 1) and charge == 4) or name == 'PP':
                                plt.plot(x,y,label = label, color = cm(neutron / 4),linewidth = 0.5)
                                ploted = True
                                plt.xlim(60, 80)
                            if (name == 'SRP2' and charge == 2) or (name =='SRP3' and charge == 3):
                                plt.plot(x,y,label = label, color = cm(neutron / 4),linewidth = 0.5)
                                ploted = True
                                plt.xlim(80, 100)
                            if (name == 'SRP6' and charge == 4) and peptide.count('o') == 1:
                                plt.plot(x,y,label = label, color = cm(neutron / 4),linewidth = 0.5)
                                ploted = True
                                plt.xlim(40, 60)
                                
                    if ploted:
                        plt.title(name +' ' + peptide)
                        plt.legend(loc = 0)
                        plt.tight_layout()
                        pp.savefig()
                        plt.close()
        pp.close()

        
        #for KB
        import matplotlib.pyplot as plt
#        import numpy as np
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(outfile)
        plt.rcParams["figure.figsize"] = [10.0,5]
        cm = plt.get_cmap('gist_rainbow')
        x = [ele[1] for ele in scan_info]
        for name in ['AGAP008922']:
#            name = name.id
            for peptide in dc_mz_neutron[name]:
                for charge in dc_mz_neutron[name][peptide]:
                    ploted = False
                    for neutron in range(neutron_max):
                        mz = dc_mz_neutron[name][peptide][charge][neutron][0]
                        percent = dc_mz_neutron[name][peptide][charge][neutron][1]
                        key = (name, peptide,charge,neutron)
                        if key in dc_mzscan and 'c' in peptide and peptide == 'ANVIDPPKVcCPDGQMLDHQGKcCRPIMG':
                            ploted = True
                            y = dc_mzscan[key]
                            label = ' charge:' +str(charge) + \
                            ' neutron' + str(neutron) +' mz:%.3f'%mz + ' ratio:%.3f' %percent
                            plt.plot(x,y,label = label, color = cm(neutron / 4),linewidth = 1)
                    if ploted:
                        plt.title(name +' ' + peptide)
                        plt.xlim(70,75)
                        plt.ylim(0,1000000)
                        plt.legend(loc = 0)
                        plt.tight_layout()
                        pp.savefig()
                        plt.close()
        pp.close()

#       for KB checking other cutting sites. 20160920
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(outfile)
        plt.rcParams["figure.figsize"] = [20.0,10]
        cm = plt.get_cmap('gist_rainbow')
        x = [ele[1] for ele in scan_info]
        for name in peptides:
            name = name.id
            for peptide in dc_mz_neutron[name]:
                for charge in dc_mz_neutron[name][peptide]:
                    ploted = False
                    mzs_hplc = []
                    for neutron in range(neutron_max):
                        mz = dc_mz_neutron[name][peptide][charge][neutron][0]
                        percent = dc_mz_neutron[name][peptide][charge][neutron][1]
                        key = (name, peptide,charge,neutron)
                        if key in dc_mzscan:
                            mzs_hplc.append(dc_mzscan[key])
                    if len(mzs_hplc) > 3:
                        mzlen = len(mzs_hplc[0])
                        mzs_hplc = np.asarray(mzs_hplc).T
                        toplot = False
                        toplotcount = 0
                        for tempn in range(mzlen):
                            if sum(mzs_hplc[tempn] >0) >2:
                                toplotcount += 1
                            if toplotcount >2:
                                toplot = True
                        if toplot:
                            for neutron in range(neutron_max):
                                mz = dc_mz_neutron[name][peptide][charge][neutron][0]
                                percent = dc_mz_neutron[name][peptide][charge][neutron][1]
                                key = (name, peptide,charge,neutron)
                                ploted = True
                                y = dc_mzscan[key]
                                label = ' charge:' +str(charge) + \
                                ' neutron' + str(neutron) +' mz:%.3f'%mz + ' ratio:%.3f' %percent
                                plt.plot(x,y,label = label, color = cm(neutron / 4),linewidth = 0.3)
                    if ploted:
                        plt.title(name +' ' + peptide)
                        plt.xlim(20,120)
                        plt.tight_layout()
                        plt.legend(loc = 0)
                        pp.savefig()
                        plt.close()
        pp.close()



    peptidesfile = r"E:\Lab\works\20160516MSMS_sytheticPeptides\analysis\10pep.fasta"
    ms1Peak = r"E:\Lab\works\20160516MSMS_sytheticPeptides\1437_XC_peptides_1.ms1.Peak"
    outfile = '20160912syntheticPeptidesMassSpec.pdf'
    dcaaModi = {'C':'d'}
    random = False
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)

    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\7pep.txt"
    ms1Peak = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.ms1.Peak"
    outfile = '20160908BarHeomolymph10kd_xlimAjusted.pdf'
    dcaaModi = {'C':'d','M':'o'}
    random = True
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)
    
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\seq4.fasta"
    ms1Peak = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.ms1.clear.Peak"
    outfile = '20160909AgKB_pep.pdf'
    dcaaModi = {'C':'c','M':'o'}
    random = True
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)
    
    #20160920 check other possible cutting sites.
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\seq_differentCutting.txt"
    from Bio import SeqIO
    peptides = list(SeqIO.parse(open(peptidesfile),'fasta'))
    fout = open(peptidesfile+'AllCutting','w')
    for peptide in peptides:
        l = len(peptide.seq)
        for n in range(0, l-18):
            seq = str(peptide.seq)[n:]
            name = peptide.id +'_' +str(len(seq))
            if len(seq) < 35:#only keep cytokine less than 35aa
                fout.write('>'+name+'\n'+seq+'\n')
    fout.close()
    
    peptidesfile = peptidesfile+'AllCutting'
    ms1Peak = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.ms1.clear.Peak"
    outfile = '20160920AgKB_pep.pdf'
    dcaaModi = {'C':'c','M':'o'}
    random = True
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)
    
    #20160922 check other possible cutting sites for Manduca SRPs
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\10pepPlus.txt"
    ms1Peak = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.ms1.Peak"
    outfile = '20160922BarHeomolymph10kd.pdf'
    dcaaModi = {'C':'d','M':'o'}
    random = True
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)
    
    #20160922 check IH
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\7pep.txt"
    ms1Peak = r"E:\store_for_D\ManducaPeptidesProteome\Raw\1380_Induced_m100.ms1.Peak"
    outfile = '20160922InducedHeomolymph_SRPs.pdf'
    dcaaModi = {'C':'c','M':'o'}
    random = True
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)
        

def mgf_dealwith20160919():
    
    import MGFparser
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\KB_Ag\1443_fullrange.mgf"
    
    filename = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\1449_10kdpep_1_200.mgf"
    targetmzcharges = [(703.821, 4),(1375.679, 2),(1045.153,3), (724.867, 4), (693.824, 4), (777.351, 4), (806.376, 3)]
    MGFparser.mzfilter(filename, targetmzcharges)
    
def cleanMS1toPeaks20160920():
    import ms1Class
    filename = r"E:\store_for_D\ManducaPeptidesProteome\Raw\1380_Control_m100.ms1"
    ms1Class.ms1ToPeaks(filename,None,10000)
    
    import ms1Class
    filename = r"E:\store_for_D\ManducaPeptidesProteome\Raw\1380_Induced_m100.ms1"
    ms1Class.ms1ToPeaks(filename,None,10000)

def cleanMS1toPeaks20161018():
    import ms1Class
    import glob
    folder = 'F:\\Insects\\2016ManducaSextaBarHemolymph10kdMs\\20161018ReRun\\ms1\\'
    filenames = glob.glob(folder + '*.ms1')
    for filename in filenames:
        ms1Class.ms1ToPeaks(filename,None,10000)


def ms1PeaksFindingSRPs20161109():
    """
    Dr. Steve Re-run my sample with fusion-ms/ms. Check these data for possible SRPs
    """
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\7pep.txt"
    ms1Peak = r"F:\Insects\2016ManducaSextaBarHemolymph10kdMs\20161018ReRun\ms1\1449rconc_alk_HH.ms1.Peak"
    outfile = '20161109BarHemolymph_Conc_HH_d.pdf'
    dcaaModi = {'C':'d','M':'o'}
    random = False
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)
    
    peptidesfile = r"E:\Lab\works\20160816BarHemolymph10kdPepNativeIdentification\7pep.txt"
    ms1Peak = r"F:\Insects\2016ManducaSextaBarHemolymph10kdMs\20161018ReRun\ms1\1449rconc_alk_HL.ms1.Peak"
    outfile = '20161109BarHemolymph_Conc_HL_d.pdf'
    dcaaModi = {'C':'d','M':'o'}
    random = False
    ms1targetFinder(peptidesfile, ms1Peak, outfile, dcaaModi,random)


def findSRPpeptidesNewWithMaxQuant20170719():
    '''
    the Bar hemolymph were re-runned. check the data with maxquant
    '''
    from pyteomics import mass
    #try to find SRP6
    seqSRP6 = 'MADSNAIVFPDEVEAAKEAAKAEVKVENGDSNTVEEAPAEPAQPLPEVALRNMIVVPPNCPPGQQMGSDGVCRVVFN'
    target = seqSRP6
    target_ms = mass.calcu
    
    file_Allpeptides = r"D:\Insects\ManducaSexta\2016ManducaSextaBarHemolymph10kdMs\20161018ReRun\raw\combined\txt\allPeptides.txt"
    import pandas as pd
    df = pd.read_table(file_Allpeptides)
    