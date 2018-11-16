#!/usr/bin/env python

import MySQLdb, subprocess, os, time, sys, multiprocessing

def connect():
    conn = MySQLdb.connect(
            host="salmiak",
            user="sparsifier",
            passwd="sparsifier",
            db="sparsifier")
    cursor = conn.cursor()
    return (conn, cursor)

def get_job():
    conn, cursor = connect()
    cursor.execute("LOCK TABLES joblist WRITE")
    cursor.execute("SELECT * FROM joblist WHERE status=0 LIMIT 1")
    rows = cursor.fetchall()

    if len(rows) > 0:
        ret = rows[0][0:2]
        cursor.execute("UPDATE joblist SET status=1 WHERE id=%d" % ret[0])
        conn.commit()
    else:
        ret = (-1, None)
        
    cursor.execute("UNLOCK TABLES")
    conn.commit()
    cursor.close()
    conn.close()

    return ret

def set_job_done(id_):
    conn, cursor = connect()
    cursor.execute("LOCK TABLES joblist WRITE")
    cursor.execute("UPDATE joblist SET status=2 WHERE id=%d" % id_)
    conn.commit()
    cursor.execute("UNLOCK TABLES")
    conn.commit()
    cursor.close()
    conn.close()

def scheduler(destdir, nthreads, libpath):
    processes = []
    ids = []
    more = True
    env = os.environ
    if libpath is not None and len(libpath) > 0:
        env["LD_LIBRARY_PATH"] = libpath
    dumpname = "input-%d.txt" % os.getpid()
    devnull = open(os.devnull, 'wb');
    while more or len(processes) > 0:
        if more and len(processes) < nthreads:
            id_, job = get_job()
            if id_ == -1:
                more = False
            else:
                open(dumpname, "w").write(job)
                ids.append(id_)
                processes.append(subprocess.Popen(
                    ["./sparsifier", destdir, dumpname, "1"],
                    stdout=sys.stdout, stderr=devnull, env=env))
        i = 0
        while i < len(processes):
            if processes[i].poll() is not None:
                set_job_done(ids[i])
                del processes[i]
                del ids[i]
            i += 1

        time.sleep(1)

def main(argv):
    if len(argv) < 2:
        return -1
    destdir = argv[1]
    libpath = None
    cpus = multiprocessing.cpu_count()
    if len(argv) > 2:
        cpus = int(argv[2])
    if len(argv) > 3:
        libpath = argv[3]
    scheduler(destdir, cpus, libpath)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

