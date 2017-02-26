import numpy as np
import scipy
import os
import sys
from numpy import linalg as LA
from time import *
import shutil
import codecs, commands
import threadpool

def call_fun(a, b, savefile):
    rc, out = commands.getstatusoutput("./dtw " + str(a) + " " + str(b) + " 15")
    print out
    cmd = "echo %s >> %s" % (out, savefile)
    os.system(cmd)

def compute_sim(indir, save_dir):
    filenames = os.listdir(indir)

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    fun_paras = []
    for i in range(len(filenames)):
        a = indir + os.sep + filenames[i]
        for j in range(0, len(filenames)):
            b = indir + os.sep + filenames[j]
            save_file = save_dir + os.sep + filenames[i].split(".")[0] + "+" + filenames[j].split(".")[0] + ".csv"
            if not os.path.exists(save_file):
                fun_paras.append(([a, b, save_file], None))
    pool = threadpool.ThreadPool(30)
    requests = threadpool.makeRequests(call_fun, fun_paras)
    [pool.putRequest(req) for req in requests]
    pool.wait()

def compute_all_sim():
    indirs = ["CENS_38_csv", "CENS_80_csv", "chroma2_38_csv", "chroma2_80_csv", "chroma_38_csv",
     "chroma_80_csv", "waon_ncp_38_csv", "widi_ncp_38_csv", "widi_ncp_80_csv"]
    save_dirs = ["CENS_38_sim", "CENS_80_sim", "chroma2_38_sim", "chroma2_80_sim", "chroma_38_sim",
    "chroma_80_sim", "waon_ncp_38_sim", "widi_ncp_38_sim", "widi_ncp_80_sim"]

    for i in range(len(indirs)):
        indir, outdir = indirs[i], save_dirs[i]
        print "start to tackle with " + indir
        compute_sim(indir, outdir)

def compute_map(indir):
    fnames, fcount = dict(), 0
    for f in os.listdir(indir):
        a, b = f.split(".")[0].split("+")
        if not a in fnames.keys():
            fnames[a] = fcount
            fcount += 1
        if not b in fnames.keys():
            fnames[b] = fcount
            fcount += 1

    sim = np.zeros((fcount, fcount))
    mark = np.array([([False]*fcount) for i in range(fcount)])
    for f in os.listdir(indir):
        a, b = f.split(".")[0].split("+")
        cmd = "cat " + indir + os.sep + f
        rc, out = commands.getstatusoutput(cmd)
        sim[fnames[a], fnames[b]] = float(out)
        if a.split("_")[0] == b.split("_")[0]:
            mark[fnames[a], fnames[b]] = True

    rank = np.argsort(sim, axis=1)
    scores = []
    for i in range(fcount):
        cover_count, all_count, score = 0, 0, 0.0
        for j in range(fcount):
            v = rank[i, j]
            if i == v:
                continue
            all_count += 1
            if mark[i, v]:
                cover_count += 1
                score += cover_count*1.0/all_count
        if cover_count > 0:
            score /= cover_count
        scores.append(score)

    ava_score = np.mean(scores)
    print indir + ": " + repr(ava_score)
    return ava_score, sim

def compute_all_map():
    indirs = ["widi_ncp_80_sim", "widi_ncp_38_sim", "waon_ncp_38_sim", "chroma_80_sim", "chroma_38_sim", "chroma2_80_sim", "chroma2_38_sim", "CENS_38_sim", "CENS_80_sim"]
    for indir in indirs:
        ava_score, sim = compute_map(indir)
        np.savetxt(fname=indir+".csv", X=sim, fmt="%.6f", delimiter=',')

def compute_AP5(indir):
    fnames, fcount = dict(), 0
    for f in os.listdir(indir):
        a, b = f.split(".")[0].split("+")
        if not a in fnames.keys():
            fnames[a] = fcount
            fcount += 1
        if not b in fnames.keys():
            fnames[b] = fcount
            fcount += 1

    sim = np.zeros((fcount, fcount))
    mark = np.array([([False]*fcount) for i in range(fcount)])
    for f in os.listdir(indir):
        a, b = f.split(".")[0].split("+")
        cmd = "cat " + indir + os.sep + f
        rc, out = commands.getstatusoutput(cmd)
        sim[fnames[a], fnames[b]] = float(out)
        if a.split("_")[0] == b.split("_")[0]:
            mark[fnames[a], fnames[b]] = True

    rank = np.argsort(sim, axis=1)
    scores = []
    cnt = 0
    for i in range(fcount):
        cover_count, all_count, score = 0, 0, 0.0
        for j in range(5):
            v = rank[i, j]
            if i == v:
                continue
            all_count += 1
            if mark[i, v]:
                cover_count += 1
                score += cover_count*1.0/all_count
        if cover_count > 0:
            score /= cover_count
            cnt += 1
        scores.append(score)

    ava_score = np.mean(scores)
    print indir + ": " + repr(ava_score) + " " + repr(cnt)
    return ava_score, sim

def compute_all_ap5():
    indirs = ["chroma_38_sim", "chroma2_80_sim", "chroma2_38_sim", "CENS_38_sim", "CENS_80_sim"]
    for indir in indirs:
        ava_score, sim = compute_AP5(indir)

def usage():
    if sys.argv[1] == "-compute_all_sim":
        compute_all_sim()
    elif sys.argv[1] == "-compute_all_map":
        compute_all_map()
    elif sys.argv[1] == "-compute_all_ap5":
        compute_all_ap5()

if __name__ == "__main__":
    usage()
