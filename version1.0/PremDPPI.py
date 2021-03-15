#!/usr/bin/python
# coding=utf-8
from string import ascii_uppercase
from string import ascii_lowercase
from collections import defaultdict,Counter
from glob import glob
import datetime,time
import pandas as pd
import os, sys, os.path, re, getopt 
import rpy2.robjects as robjects

ascii_cases = ascii_uppercase + ascii_lowercase  

r = robjects.r

normal_format_pro = ['CYS','GLN','ILE','SER','VAL','MET','ASN','PRO','LYS','THR','PHE','ALA','HIS','GLY','ASP','LEU','ARG','TRP','GLU','TYR']

# map residue name three letters to one
map_three_to_one = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
       "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
       "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
       "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
       "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
       "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
       "PTR": "X", "XLE": "X", "XAA": "X" }

# Amino acid masses. http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
mapm = {'A':71.0788,'R':156.1875,'N':114.1038,'D':115.0886,'C':103.1388,'E':129.1155,'Q':128.1307,'G':57.0519,'H':137.1411,'I':113.1594,'L':113.1594,'K':128.1741,'M':131.1926,'F':147.1766,'P':97.1167,'S':87.0782,'T':101.1051,'W':186.2132,'Y':163.176,'V':99.1326}

# set up path.
workdir = Your working directory
pathvmd = path for running VMD software  # /usr/local/bin/vmd
pathcharmm = path for running CHARMM software  # /usr/local/bin/charmm
pathnamd2 = path for running NAMD software  # /usr/local/bin/NAMD_2.12_Source/Linux-x86_64-g++/namd2
pathmkdssp = path for running DSSP software  # /usr/local/bin/mkdssp
pathpsiblast = path for running PSI-BLAST software  # /usr/local/bin/blast/psiblast
pathblastdb = path for blastdb  # /usr/local/bin/blastdb/nr
pathprovean = path for PROVEAN software  # /usr/bin/provean.sh
pathipot = path for iPot software  # /usr/local/bin/iPot/

pathpara = workdir + "inputfiles"

# input jobid
jobid = ''
myopts, args = getopt.getopt(sys.argv[1:], "i:z")
for o, a in myopts:
    if o == '-i':
        jobid = a
        # print jobid
    else:
        print("Usage: %s -i jobid" % sys.argv[0])

# set up input and output file name for each job
jobpath = workdir + jobid
pathoutput = workdir + jobid + "out"

os.system("mkdir %s" % pathoutput)
in_file = jobpath + '/' + jobid + '.input'  # inputfile name
out_file = jobpath + '/' + jobid + '.classification'  # outputfile name


# process biounit model and write another colummn for indicating chains information from inputfiles
def ProPDB1():
    pdball = []
    with open(in_file, 'r') as f:
        next(f)
        for line in f:
            ff = line.split('\t')
            pdb = ff[0]  # 1AUD.pdb or 1AUD.pdb1
            usedchains = ff[1]+ff[2]
            if pdb not in pdball:
                pdball.append(pdb)
                if pdb.split(".")[1] != "pdb":
                    ST_play = False
                    ffpdb = open(pathoutput + '/' + pdb.split(".")[0].lower() + '_p.pdb', 'w')
                    fpdb = open(jobpath + '/' + pdb, 'r')
                    for linepdb in fpdb:
                        if linepdb[0:5] == "MODEL":
                            CountModel = linepdb.split()[1]
                            ST_play = True
                            continue
                        if ST_play:
                            if linepdb[:4] == "ATOM" and linepdb[21] in usedchains and linepdb[17:20].strip() in normal_format_pro: 
                                ffpdb.write("%s                  %s\n" % (
                                linepdb[0:54].strip('\r\n'), str(linepdb[21:22]) + '_' + str(CountModel)))
                else:
                    ffpdb = open(pathoutput + '/' + pdb.split(".")[0].lower() + '_p.pdb', 'w')  # 1yag_p.pdb
                    fpdb = open(jobpath + '/' + pdb, 'r')
                    ST_play = True
                    for linepdb in fpdb:
                        line_list = re.split(r'\s+', linepdb)
                        if (line_list[0] == 'MODEL') and ST_play:
                            countmodel = line_list[1]
                            ST_play = False
                        if (line_list[0] == 'MODEL') and (line_list[1] != countmodel):
                            break
                        if linepdb[:4] == 'ATOM' and linepdb[21] in usedchains and linepdb[17:20].strip() in normal_format_pro:  
                            ffpdb.write("%s                  %s\n" % (linepdb[0:54].strip('\r\n'), str(linepdb[21:22]) + '_' + str(1)))
            else:
                continue
    ffpdb.close()
    fpdb.close()


# produce 1yag_p.pdb
def del_unknown_incomplete():
    residues = [k for k, v in map_three_to_one.items() if v!= 'X']
    f = open(in_file, 'r')
    next(f)
    for line in f:
        ff = line.split("\t")
        pdbid = ff[0].split('.')[0].lower()
        pdbfile = open('{}/{}_p.pdb'.format(pathoutput,pdbid)).readlines()
        # delete unknow resdues
        pdbfile_del_unknown = [i for i in pdbfile if i[17:20] in residues]
        # delete incomplete residues
        final_row = pdbfile_del_unknown[-1]
        last = ''
        above = []
        allresidues = []
        for row in pdbfile_del_unknown:     
            if row[17:26] == last and row == final_row:
                above.append(row)
                atoms = [i[13:16].strip() for i in above]
                if set(['C','N','O','CA']).issubset(set(atoms)):
                    allresidues.append(above)
            elif row[17:26] == last and row != final_row:
                above.append(row)
            else:                            # when read different residue
                if len(above)>=4:
                    atoms = [i[13:16].strip() for i in above]
                    if set(['C','N','O','CA']).issubset(set(atoms)):
                        allresidues.append(above)
                above = [row] 
            last = row[17:26]
        # write out
        with open('{}/{}_p_test.pdb'.format(pathoutput,pdbid),'w') as fw:
            fw.write(''.join([y for x in allresidues for y in x]))
        break
    os.system('mv {}/{}_p_test.pdb {}/{}_p.pdb'.format(pathoutput,pdbid, pathoutput,pdbid))
    f.close()


# split chains and produce pdb files for each chain of wild type pdb.
def splitchain():
    pdball = []
    f = open(in_file, 'r')
    next(f)
    for line in f:
        ff = line.split("\t")
        pdbfile = ff[0]
        pdb = pdbfile.split(".")[0].upper()
        partner1 = ff[1].split(".")
        partner2 = ff[2].split(".")
        if pdbfile not in pdball:
            for chains in list(partner1 + partner2):
                f1 = open(pathoutput + '/' + pdb.lower() + '_p.pdb', 'r')
                fw = open(pathoutput + '/' + pdb + '_' + chains + '.pdb', 'w')
                for line1 in f1:
                    if chains == line1[72:].strip():
                        fw.write(line1)
                f1.close()
                fw.close()
            pdball.append(pdbfile)
        else:
            continue
    f.close()


# processing pdb files and prepair input files for provean
def CleanPdb():
    first_line = open(in_file, 'r').readlines()[0][:-1]
    fw = open(in_file + ".cleaned", "w")
    fw.write("%s\t%s\t%s\t%s\t%s\n" % (first_line, "PDBid", "NewPartner1", "NewPartner2", "Mutation_cleaned"))

    second_line = open(in_file, 'r').readlines()[1]
    ff = second_line.split("\t")
    pdb = ff[0].split(".")[0].upper()
    partner1 = ff[1].split(".")
    partner2 = ff[2].split(".")

    mapchainarray = []
    counti = 0
    for chains in list(partner1 + partner2):
        cc = (chains, ascii_cases[counti])
        mapchainarray.append(cc)
        counti += 1
        mapchaindict = dict(iter(mapchainarray))

    newpartner1 = ''
    for chains in list(partner1):
        newpartner1 += mapchaindict[chains]
    newpartner2 = ''
    for chains in list(partner2):
        newpartner2 += mapchaindict[chains]

    countchain = 1
    for chains in list(partner1 + partner2):
        fwpdb = open(pathoutput + "/" + pdb + "_CH" + str(countchain) + ".pdb", "w")
        fvar = open(pathoutput + "/" + pdb + "_" + chains + ".var", "w")
        countchain += 1
        count = 1
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        resname = fpdb.readlines()[0][17:20].strip()
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        resnum = fpdb.readlines()[0][22:27].strip()
        line1 = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r").readlines()[0]
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        atomname = fpdb.readlines()[0][13:16].strip()
        fwpdb.write("%s %s%s   %s %s" % (line1[0:16], line1[17:21], mapchaindict[chains], str(count), line1[27:]))

        mutchainall = []
        f = open(in_file, 'r')
        _unused = next(f)
        for line in f:
            ff = line.split("\t")
            mut = ff[4].upper()
            mutchain = ff[3]
            if str(chains) == str(mutchain) and mutchain not in mutchainall:
                fseq = open(pathoutput + "/" + pdb + "_" + mutchain + ".seq", "w")
                fseq.write("%s %s\n" % (">", pdb + mutchain))
                fseq.write("%s" % (map_three_to_one[resname]))
            mutchainall.append(mutchain)
            if str(resnum) == str(mut[1:-1]) and str(map_three_to_one[resname]) == str(mut[0:1]) and str(chains) == str(mutchain):
                fw.write("%s\t%s\t%s\t%s\t%s\n" % (line[:-1], pdb, newpartner1, newpartner2, str(map_three_to_one[resname]) + mapchaindict[chains] + str(count) + mut[-1:]))
                fvar.write("%s\n" % (str(map_three_to_one[resname]) + str(count) + mut[-1:]))
        f.close()

        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        for linepdb in fpdb.readlines()[1:]:
            resnamepdb = linepdb[17:20].strip()
            resnumpdb = linepdb[22:27].strip()
            atomnamepdb = linepdb[13:16].strip()
            if resnamepdb == resname and resnumpdb == resnum:
                if atomnamepdb != atomname:
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 10 and count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 100 and count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 1000 and count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                else:
                    continue
            else:
                if atomnamepdb != atomname:
                    count += 1
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 10 and count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 100 and count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 1000 and count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))

                    mutchainall = []
                    with open(in_file, 'r') as f:
                        _unused = next(f)
                        for line in f:
                            ff = line.split("\t")
                            mut = ff[4]
                            mutchain = ff[3]
                            if str(chains) == str(mutchain) and mutchain not in mutchainall:
                                fseq.write("%s" % (map_three_to_one[resnamepdb]))
                            mutchainall.append(mutchain)
                            if str(resnumpdb) == str(mut[1:-1]) and str(map_three_to_one[resnamepdb]) == str(mut[0:1]) and str(chains) == str(mutchain):
                                fw.write("%s\t%s\t%s\t%s\t%s\n" % (
                                line[:-1], pdb, newpartner1, newpartner2, str(map_three_to_one[resnamepdb]) + mapchaindict[chains] + str(count) + mut[-1:]))
                                fvar.write("%s\n" % (str(map_three_to_one[resnamepdb]) + str(count) + mut[-1:]))
                else:
                    continue
            resname = linepdb[17:20].strip()
            resnum = linepdb[22:27].strip()
            atomname = linepdb[13:16].strip()
        fpdb.close()
        fwpdb.close()
        fvar.close()
    fseq.close()
    fw.close()


# produce wild type pdb structure for DSSP running.
def wtpdb():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        if pdb not in pdball:
            os.system('cat %s/%s_CH*.pdb > %s/%s.pdb' % (pathoutput, pdb, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


# produce psf and pdb files of wild-type with vmd.
def vmd_wt():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    template = open(pathpara + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = ff[8]
        partner2 = ff[9]
        NumChain = int(len(partner1 + partner2))
        if pdb not in pdball:
            vmd_pdb = template.replace('protname', pdb).replace('NumChain', str(NumChain)).replace('pathinput', pathpara).replace('pathoutput', pathoutput)
            with open(pathoutput + '/vmd_' + pdb + '.pgn', 'w') as fw:
                fw.write(vmd_pdb)
            os.system('%s -dispdev text -e %s/vmd_%s.pgn' % (pathvmd, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


# produce charmm psf and crd files for individual chains
def charmmfile_wt():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        NumP = int(len(partner1 + partner2))
        countp = 1
        if pdb not in pdball:
            countp = 1
            for chains in (list(partner1) + list(partner2)):
                os.system('grep "^.\{21\}%s" %s/%s_vmd.pdb > %s/%s_%s.pdb' % (chains, pathoutput, pdb.upper(), pathoutput, pdb, 'ch' + str(countp)))
                ffpdb = open(pathoutput + '/' + pdb + '_ch' + str(countp) + '.pdb', 'a')
                ffpdb.write("%s" % ('END'))
                ffpdb.close()
                os.system('%s <%s/setup.inp path=%s protname=%s I=%s pathpara=%s> %s/setup_%s.out' % (pathcharmm, pathpara , pathoutput, pdb, countp, pathpara, pathoutput, pdb + '_ch' + str(countp)))
                countp += 1
            os.system('%s <%s/append.inp path=%s protname=%s NumP=%s pathpara=%s> %s/append_%s.out' % (pathcharmm, pathpara, pathoutput, pdb, NumP, pathpara, pathoutput, pdb))
            os.system('%s <%s/charmm_to_namd.inp path=%s protname=%s pathpara=%s> %s/charmm_to_namd_%s.out' % (pathcharmm, pathpara, pathoutput, pdb, pathpara, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


# get interface for wild type. interface control.
def interface():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        partner1 = ff[8]
        partner2 = ff[9]
        NumP1 = int(len(partner1))
        NumP2 = int(len(partner2))
        if pdb not in pdball:
            p1 = 'ch1'
            p2 = 'ch' + str(NumP1 + 1)
            for countp1 in range(2, NumP1 + 1):
                p1 += ' .or. segid ch' + str(countp1)
            for countp2 in range(NumP1 + 2, NumP1 + NumP2 + 1):
                p2 += ' .or. segid ch' + str(countp2)
            os.system('sed -e "s/p1/%s/g" -e "s/p2/%s/g" %s/interface.inp > %s/interface_%s.inp' % (p1, p2, pathpara, pathoutput, pdb))
            os.system('%s <%s/interface_%s.inp path=%s protname=%s pathpara=%s> %s/interface_%s.out' % (pathcharmm, pathoutput, pdb, pathoutput, pdb, pathpara, pathoutput, pdb))
            os.system('grep \'DELTAMS <\' %s/interface_%s.out | cut -c 24-100 > %s/dsasa_%s.out' % (pathoutput, pdb, pathoutput, pdb))
            fdsasa = open(pathoutput + '/dsasa_' + pdb + '.out', 'r')
            global dsasa
            dsasa = float(fdsasa.readline()[1:-2])
            pdball.append(pdb)
    f.close()


# run rest.png (change directory) for preparing the restrain pdb file for minimization of NAMD running.
def runrest():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    template = open(pathpara + '/rest.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        if pdb not in pdball:
            rest_wild = template.replace('protname', pdb).replace('pathoutput', pathoutput)
            with open(pathoutput + '/rest_' + pdb + '.pgn', 'w') as fw:
                fw.write(rest_wild)
            os.system('%s -dispdev text -e %s/rest_%s.pgn' % (pathvmd, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


# run minimization.
def minimization():
    pdball = []
    runstep = "100"
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    template = open(pathpara + '/min.conf').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        if pdb not in pdball:
            min_wild = template.replace('protname', pdb).replace('pathoutput', pathoutput).replace('runstep', runstep).replace('pathinput', pathpara)
            with open(pathoutput + '/min_' + pdb + '.conf', 'w') as fw:
                fw.write(min_wild)
            os.system('%s %s/min_%s.conf' % (pathnamd2, pathoutput, pdb))
        else:
            continue
    f.close()


def frames():
    """
    get frame from minimization. 
    """
    pdball = []
    startframe = 100
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        if pdb not in pdball:
            os.system('%s <%s/frame.inp path=%s protname=%s startframe=%s pathinput=%s> %s/frame_%s.out' % (pathcharmm, pathpara, pathoutput, pdb, startframe, pathpara, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


def get_p2_len():
    f = open(in_file + ".cleaned", 'r')
    next(f)
    fw = open(pathoutput+'/'+jobid+'_len_p2.txt','w')
    fw.write('PDBid\tMutation_cleaned\tL_p2\n')
    for line in f:
        ff = line.strip().split("\t")
        pdb = ff[7].lower()  # 1yag
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10]
        mutchain = mut[1]  # a

        if mutchain.upper() in partner1:
            partner = 'p1'
            op = partner2
        else:
            partner = 'p2'
            op = partner1

        count_op = set()
        with open(pathoutput+ '/' + pdb.upper()+'.pdb') as fpdb:
            for ffpdb in fpdb:
                chain = ffpdb[21]
                if chain in op:
                    count_op.add((ffpdb[21], ffpdb[22:27]))
        len_op = len(count_op)
        fw.write('%s\t%s\t%s\n'%(pdb.upper(),mut,len_op))


def RunProvean():
    pdbchain = []
    f = open (in_file+".cleaned",'r')
    next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        chain = ff[3]
        if pdb+chain not in pdbchain:
            os.system('%s -q %s/%s_%s.seq -v %s/%s_%s.var > %s/provean_%s_%s.out --num_threads 30' % (pathprovean, pathoutput,pdb,chain,pathoutput,pdb,chain,pathoutput,pdb,chain))
            pdbchain.append(pdb+chain)
    f.close()


def get_provean_score():
    f = open(in_file + ".cleaned", 'r')
    next(f)
    fw = open(pathoutput+'/'+jobid+'_provean_score.txt','w')
    first_line = open(in_file + ".cleaned", 'r').readlines()[0][:-1]
    fw.write('%s\t%s\n'%(first_line,'DCS')) 
    for line in f:
        ff = line.strip().split("\t")
        pdb = ff[7].lower()  # 1yag
        mut = ff[10].lower()  # da11a
        ST_play = False
        fp = open (pathoutput+"/provean_"+pdb.upper()+"_"+ff[3]+".out","r")
        for provean in fp:
            ffp = provean.split("\t")
            if provean[2:11] == "VARIATION":
                ST_play = True
                continue
            if ST_play:
                if ffp[0] == (mut[0] +mut[2:]).upper():
                    provean_score = ffp[1][:-1]
                    fw.write('%s\t%s\n'%(line[:-1],provean_score))
    fw.close()


def dssp():
    """ produce dssp files showing secondary structure infromation for crystal structure of wild-type PDB. """
    pdball = []
    mutchainall = []
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        chain = ff[10][1:2]
        mutchain = ff[3]  # A_1
        if pdb not in pdball:
            # COMPLEX
            os.system('%s -i %s/%s.pdb -o %s/%s.dssp' % (pathmkdssp, pathoutput, pdb, pathoutput, pdb))
            pdball.append(pdb)

    f.close()


# calculate sasa for each residue
def surfacefile():
    f = open(in_file + ".cleaned", 'r')
    next(f)
    template = open(pathpara + '/surface.inp').read()
    for line in f:
        ff = line.split("\t")
        mutchain = ff[3].lower()
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        chain = ff[10][1:2]
        resnum = ff[10][2:-2]
        partner1 = ff[8]
        partner2 = ff[9]
        NumP1 = int(len(partner1))
        NumP2 = int(len(partner2))
        p1 = 'ch1'
        p2 = 'ch' + str(NumP1 + 1)
        for countp1 in range(2, NumP1 + 1):
            p1 += ' -\n .or. segid ch' + str(countp1)
        for countp2 in range(NumP1 + 2, NumP1 + NumP2 + 1):
            p2 += ' -\n .or. segid ch' + str(countp2)

        # wild
        fwt = open(pathoutput + "/surface_" + protname + "_wt.inp", "w")
        fwt_surface = template.replace('p1', p1).replace('p2', p2)
        fwt.write(fwt_surface)
        fwt.close()
    f.close()


def surface():
    f = open(in_file + ".cleaned", 'r')
    next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        chain = ff[10][1:2]
        resnum = ff[10][2:-2]
        partner1 = ff[8]
        partner2 = ff[9]
        NumP1 = int(len(partner1))
        NumP2 = int(len(partner2))
        ST_play = False
        mapchainarray = []
        count = 1
        for chains in (list(partner1) + list(partner2)):
            cc = (chains, 'ch' + str(count))
            mapchainarray.append(cc)
            count += 1
        mapchaindict = dict(iter(mapchainarray))
        # wild
        if chain in list(partner1):
            os.system('%s <%s/surface_%s.inp pathinput=%s path=%s protname=%s chainm=%s resnumber=%s partner=%s> %s/surface_%s.out' % (pathcharmm, pathoutput, protname+'_wt', pathpara, pathoutput, pdb, mapchaindict[chain], resnum, "a", pathoutput, protname+'_wt'))
        else:
            os.system('%s <%s/surface_%s.inp pathinput=%s path=%s protname=%s chainm=%s resnumber=%s partner=%s> %s/surface_%s.out' % (pathcharmm, pathoutput, protname+'_wt', pathpara, pathoutput, pdb, mapchaindict[chain], resnum, "b", pathoutput, protname+'_wt'))

    f.close()


def judge_interface_residue():
    f = open(in_file + ".cleaned", 'r')
    next(f)
    fx = open(pathoutput+'/'+jobid+'_interface.txt','w')
    first_line = open(in_file + ".cleaned", 'r').readlines()[0][:-1]
    fx.write('%s\t%s\n'%(first_line,'Interface'))
    for line in f:
        ff = line.strip().split("\t")
        pdb = ff[7].lower()  # 1yag
        Mutation_PDB = ff[4]  # D11A
        mut = ff[10].lower()  # da11a
        protname = pdb + '_' + mut  # 1yag_da11a

        # add location
        ## wild
        rasac_wild = float(os.popen('grep \' RASAC <\' %s/surface_%s.out | cut -c 22-100' % (pathoutput, protname+'_wt')).read().strip()[1:-1])
        rasamp_wild = float(os.popen('grep \' RASAMP <\' %s/surface_%s.out | cut -c 23-100' % (pathoutput, protname+'_wt')).read().strip()[1:-1])
        drasap_wild = float(os.popen('grep \' DRASAP <\' %s/surface_%s.out | cut -c 23-100' % (pathoutput, protname+'_wt')).read().strip()[1:-1])
        ## if mutation on interface
        locationx = ''
        if drasap_wild > 0.0 and rasamp_wild > 0.25 and rasac_wild < 0.25:
            locationx = "yes"
        if drasap_wild > 0.0 and rasac_wild > 0.25:
            locationx = "yes"
        if drasap_wild > 0.0 and rasamp_wild < 0.25:
            locationx = "yes"
        if drasap_wild == 0.0 and rasac_wild > 0.25:
            locationx = "no"
        if drasap_wild == 0.0 and rasac_wild < 0.25:
            locationx = "no"
        fx.write('%s\t%s\n'%(line[:-1],locationx))

    fx.close()
    f.close()


# Amino acid composition-mutchain  , Solart: A Structure-Based Method To Predict Protein Solubility And Aggregation. 2019.BioRxiv 
def AAcompo():
    residues = ['T', 'I', 'V', 'P', 'A', 'L', 'D', 'E', 'R', 'K', 'F', 'W', 'Y']
    namelist_wt = ['AAcompo_'+i+'_wt' for i in residues+['RKDE','FWY','IVLAP']]

    with open(in_file+'.cleaned') as f, open(pathoutput+'/amino_acid_composition_mutchain.txt','w') as fw:
        fw.write('PDBid\tPartner1\tPartner2\tMutChain\tMutation_PDB\tMutation_cleaned\t'+'\t'.join(namelist_wt)+'\n')
        next(f)
        pdball = []
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[7]
            partner1 = linelist[1]
            partner2 = linelist[2]
            mutchain = linelist[3]
            mutation_pdb = linelist[4]
            mut = linelist[-1]
            location = int(mut[2:-1])
            # wild
            if pdb+mutchain not in pdball:
                pdball.append(pdb+mutchain)
                seqwt = list(open('{}{}out/{}_{}.seq'.format(workdir,jobid,pdb,mutchain)).readlines()[1])
                length = len(seqwt)
                composition_wt = dict(Counter(seqwt))
                dwt = {'AAcompo_'+i+'_wt':((composition_wt[i]+0.00)/length) if i in composition_wt.keys() else 0.00 for i in residues}
                dwt.update({'AAcompo_RKDE_wt':dwt['AAcompo_R_wt']+dwt['AAcompo_K_wt']+dwt['AAcompo_D_wt']+dwt['AAcompo_E_wt']})
                dwt.update({'AAcompo_FWY_wt':dwt['AAcompo_F_wt']+dwt['AAcompo_W_wt']+dwt['AAcompo_Y_wt']})
                dwt.update({'AAcompo_IVLAP_wt':dwt['AAcompo_I_wt']+dwt['AAcompo_V_wt']+dwt['AAcompo_L_wt']+dwt['AAcompo_A_wt']+dwt['AAcompo_P_wt']})

            fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pdb, partner1, partner2, mutchain, mutation_pdb, mut, '\t'.join([str(dwt[i]) for i in namelist_wt])))


def getdsspRatio(dsspfile):
    flag = False
    total_res = set()
    d = defaultdict(set)
    result = defaultdict(float)
    with open(dsspfile) as fdssp:
        for ldssp in fdssp:
            if ldssp.startswith('  #  RESIDUE AA STRUCTURE'):
                flag = True
            if flag:
                resnum = ldssp[5:10].strip()
                chain = ldssp[11]
                ss = ldssp[16]
                total_res.add(resnum+chain)

                # α-helix
                if ss == 'H':
                    d['H'].add(resnum+chain)
                # result
                result['H'] = float(len(d['H']))/len(total_res)
    return result


def get_secondary_proportion():
    # get result for each mutation
    d = defaultdict(lambda:defaultdict(lambda:defaultdict(float)))
    with open(in_file+'.cleaned') as f:
        next(f)
        for line in f:
            pdb = line.strip().split('\t')[0].split('.')[0]
            mut = line.strip().split('\t')[-1]
            # wild
            d[pdb]['wild']=getdsspRatio(pathoutput+'/'+pdb+'.dssp')

    # write to outfile
    with open(pathoutput+'/'+jobid+'_Secondary_Structure_Proportion.txt','w') as fw:
        fw.write('PDBid\tP_helix\n')
        for k1,v1 in d.items():
            fw.write(k1+'\t'+str(d[k1]['wild']['H'])+'\n')


def iPot():
    potential_list = ['rrce20']
    f = open(in_file + '.cleaned').readlines()[1:]
    pdb = f[0].split('\t')[0].split('.')[0].lower()  # 1YAG
    # wile type
    wt_res = []
    for i in potential_list:
        os.system('%s -r %s > %s' % (pathipot+i, pathoutput+'/'+pdb+'_min.pdb',pathoutput+'/'+pdb+'_min_'+i+'.log'))
        wt_res.append(open(pathoutput+'/'+pdb+'_min_'+i+'.log').readlines()[0].strip().split('=')[-1].strip())
    with open(pathoutput+'/'+jobid+'_iPot.txt','w') as fw:
        fw.write('\t'.join(['PDBid','E_concat'])+'\n')
        fw.write('\t'.join([pdb.upper()]+wt_res)+'\n')



def get_model_features():
    provean_f = pd.read_csv(pathoutput+'/'+jobid+'_provean_score.txt', sep='\t',dtype=str)
    aacom_f = pd.read_csv(pathoutput+'/amino_acid_composition_mutchain.txt', sep='\t',dtype=str)
    aacom_f.rename(columns={'AAcompo_RKDE_wt':'P_RKDE','AAcompo_FWY_wt':'P_FWY','AAcompo_IVLAP_wt':'P_IVLAP','AAcompo_T_wt':'P_T'},inplace=True)
    out_f1 = pd.merge(provean_f,aacom_f,on=['PDBid','Partner1','Partner2','MutChain','Mutation_PDB','Mutation_cleaned'],how='left')
    helix_f = pd.read_csv(pathoutput+'/'+jobid+'_Secondary_Structure_Proportion.txt', sep='\t',dtype=str) 
    out_f2 = pd.merge(out_f1,helix_f,on=['PDBid'],how='left')
    hydro_f = pd.read_csv(pathpara+'/hydrophobicity_scales.txt',sep='\t',dtype=str)
    hydro_f['residue'] = [map_three_to_one[i] for i in hydro_f['Residue']]
    dict_hydro = hydro_f.groupby('residue')['OMH'].apply(list).to_dict()
    out_f2['DOMH'] = [float(','.join(dict_hydro[x[0]]))-float(','.join(dict_hydro[x[-1]])) for x in out_f2['Mutation_PDB']]
    out_f2['DMass'] = [mapm[x[0]]-mapm[x[-1]] for x in out_f2['Mutation_PDB']]
    p2_f = pd.read_csv(pathoutput+'/'+jobid+'_len_p2.txt', sep='\t',dtype=str)
    out_f3 = pd.merge(out_f2,p2_f,on=['PDBid','Mutation_cleaned'],how='left')
    ipot_f = pd.read_csv(pathoutput+'/'+jobid+'_iPot.txt', sep='\t',dtype=str)
    interface_f = pd.read_csv(pathoutput+'/'+jobid+'_interface.txt', sep='\t',dtype=str)
    out_f4 = pd.merge(out_f3,interface_f,on=['PDBfile','Partner1','Partner2','MutChain','Mutation_PDB','Result_Id','isPI','PDBid','NewPartner1','NewPartner2','Mutation_cleaned'],how='left')
    merge_f = pd.merge(out_f4,ipot_f,on=['PDBid'],how='left')
    merge_f1 = merge_f[['PDBfile','Partner1','Partner2','MutChain','Mutation_PDB','Result_Id','isPI','PDBid','NewPartner1','NewPartner2','Mutation_cleaned','Interface','DCS','P_RKDE','P_FWY','P_IVLAP','P_T','P_helix','DOMH','DMass','L_p2','E_concat']].drop_duplicates()
    merge_f1.to_csv(in_file + '.cleaned.outdata',index=0,sep='\t')


def Prediction():
    outdata = in_file + '.cleaned.outdata'
    robjects.globalenv["outdata"] = outdata
    robjects.globalenv["workdir"] = workdir
    robjects.globalenv["jobpath"] = jobpath
    robjects.globalenv["jobid"] = jobid
    r('''library(randomForest)''')
    r('''test = read.table(outdata,header=T,sep="\t")''')
    r('''filename_rf = paste(workdir, 'inputfiles/PremPPI.RData',sep = '')''')
    r('''load(file = filename_rf)''')
    r('''test$PremDPPI = predict(premppi.rf,test)''')   # model
    r('''write.table(test,file = paste(jobpath,'/',jobid , '.input.output',sep=''),quote = FALSE,col.names = TRUE,row.names = FALSE,sep = '\t')''') 
    
    fout = pd.read_csv(jobpath + '/' + jobid + '.input.output',sep='\t',dtype=str)
    fout1 = fout[['PDBfile','Partner1','Partner2','MutChain','Mutation_PDB','Result_Id','isPI','Interface','PremDPPI']].drop_duplicates()
    fout1.to_csv(out_file,sep='\t',index=0)
    os.system('rm %s' % (jobpath + '/' + jobid + '.input.output'))


def inputfoldx():
    """command=BuildModel"""
    f = open(in_file + ".cleaned", 'r')
    next(f)
    for line in f:
        ff = line.split("\t")
        mut = ff[10][:-1]  # DA11A
        mutchain = mut[1:2]  # A
        partner1 = ff[8]  # A
        pdb = ff[7]  # 1YAG
        with open('individual_list_' + jobid + '_' + mut + '.txt', 'w') as fic:
            fic.write('%s' % (mut + ';'))
        with open('foldx_buildmodel_' + jobid + '_' + mut + '.txt', 'w') as fsc:
            fsc.write('command=BuildModel\npdb=%s\nmutant-file=%s\n' % (
                jobid + '_' + mut + '.pdb', 'individual_list_' + jobid + '_' + mut + '.txt'))
        os.system("cp %s/%s.pdb %s_%s.pdb" % (pathoutput, pdb, jobid, mut))
    f.close()


def runfoldx_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        mut = ff[10][:-1]
        os.system('./foldx -f foldx_buildmodel_%s_%s.txt' % (jobid, mut))
        # delete files
        os.system("rm %s/%s_%s_temp.pdb" % (pathoutput, jobid, mut))
        os.system("rm %s_%s.pdb" % (jobid, mut))
        os.system("rm WT_%s_%s_1.pdb" % (jobid, mut))
        os.system("rm individual_list_%s_%s.txt" % (jobid, mut))
        os.system("rm foldx_buildmodel_%s_%s.txt" % (jobid, mut))
        os.system("rm Average_%s_%s.fxout" % (jobid, mut))
        os.system("rm Raw_%s_%s.fxout" % (jobid, mut))
        os.system("rm PdbList_%s_%s.fxout" % (jobid, mut))
        os.system("rm %s_%s.fxout" % (jobid, mut))
        os.system("mv Dif_%s_%s.fxout %s" % (jobid, mut, pathoutput))
        os.system("mv %s_%s_1.pdb %s/%s_%s.pdb" % (jobid, mut, pathoutput, pdb, mut))
    f.close()


def splitchain_mut():
    f = open(in_file + ".cleaned", 'r')
    next(f)
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10][:-1]
        pdb = pdb + '_' + mut
        countchain = 1
        for chains in (list(partner1) + list(partner2)):
            os.system('grep "^.\{21\}%s" %s/%s.pdb > %s/%s_%s.pdb' % (chains, pathoutput, pdb, pathoutput, pdb, 'CH' + str(countchain)))
            countchain += 1
    f.close()


# produce psf and pdb files with vmd. must do this for charmm input, Changing the file of vmd.pgn's path as you need and force field.
def vmd_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = next(f)
    template = open(pathpara + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10][:-1]
        protname = pdb + '_' + mut
        NumChain = int(len(partner1 + partner2))
        vmd_pdb = template.replace('protname', protname).replace('NumChain', str(NumChain)).replace('pathinput', pathpara).replace('pathoutput', pathoutput)
        with open(pathoutput + '/vmd_' + protname + '.pgn', 'w') as fw:
            fw.write(vmd_pdb)
        os.system('%s -dispdev text -e %s/vmd_%s.pgn' % (pathvmd, pathoutput, protname))
    f.close()


def main1():
    ProPDB1()
    del_unknown_incomplete()
    splitchain()
    CleanPdb()
    wtpdb()
    vmd_wt()
    charmmfile_wt()
    interface()
    if dsasa <= -200:  # can change
        runrest()
        minimization()
        frames()
        get_p2_len()
        RunProvean()
        get_provean_score()
        dssp()
        surfacefile()
        surface()
        judge_interface_residue()
        AAcompo()
        get_secondary_proportion()
        iPot()
        get_model_features()
        Prediction()
        # run mutant
        ispi = open(in_file, 'r').readlines()[1].split("\t")[-1][:-1]
        if ispi == '1':
            inputfoldx()
            runfoldx_mut()
            splitchain_mut()
            vmd_mut()
    else:
        f = open(out_file + '.error1', 'w')
        f.write("%s" % "Your interaction parters do not contact with each other. Dsasa > -200.")


if __name__ == '__main__':
    start =time.time()
    main1()
    end = time.time()
    print('Running time: %s Seconds'%(end-start))
