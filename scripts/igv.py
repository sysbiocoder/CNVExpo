# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 10:49:35 2021

@author: sth036
"""

import socket
import os.path as op
import os
import sys
class IGV(object):
    _socket = None
    _path = None
    def __init__(self, host='127.0.0.1', port=60151, snapshot_dir='/tmp/igv'):
        self.host = host
        self.port = port
        self.commands = []
        self.connect()
        self.set_path(snapshot_dir)

    @classmethod
    def start(cls, jnlp="igv.jnlp", url="http://www.broadinstitute.org/igv/projects/current/"):
        import subprocess
        from threading import Thread
        import time

        def readit(ffrom, fto, wait):
            for line in iter(ffrom.readline, b''):
                if "Listening on port" in line:
                    wait[0] = False
                fto.write(line + '\n')
            ffrom.close()

        p = subprocess.Popen("/usr/bin/javaws -Xnosplash %s%s" % (url, jnlp),
                             shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        wait = [True]
        _tout = Thread(target=readit, args=(p.stdout, sys.stdout, wait))
        _terr = Thread(target=readit, args=(p.stderr, sys.stderr, wait))
        _tout.daemon = _terr.deamon = True
        _tout.start()
        _terr.start()
        while p.poll() is None and wait[0]:
            time.sleep(10)
            print("waiting", wait)

    def connect(self):
        if self._socket:
            self._socket.close()
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._socket.connect((self.host, self.port))

    def go(self, position):
        return self.send('goto ' + position)
    goto = go

    def genome(self, name):
        return self.send('genome ' + name)

    def load(self, url):
        return self.send('load ' + url)

    def region(self, contig, start, end):
        return self.send(' '.join(map(str, ['region', contig, start, end])))

    def sort(self, option='base'):
        """
        options is one of: base, position, strand, quality, sample, and
        readGroup.
        """
        assert option in ("base", "position", "strand", "quality", "sample",
                          "readGroup")
        return self.send('sort ' + option)

    def set_path(self, snapshot_dir):
        if snapshot_dir == self._path:
            return
        if not op.exists(snapshot_dir):
            os.makedirs(snapshot_dir)

        self.send('snapshotDirectory %s' % snapshot_dir)
        self._path = snapshot_dir

    def expand(self, track=''):
        self.send('expand %s' % track)

    def collapse(self, track=''):
        self.send('collapse %s' % track)

    def clear(self):
        self.send('clear')

    def send(self, cmd):
        # socket in Python2 oprates with strings
        if sys.version_info.major == 2:
            self._socket.send(cmd + '\n')
            return self._socket.recv(4096).rstrip('\n')
        # while socket in Python3 requires bytes
        else:
            self.commands.append(cmd)
            cmd = cmd + '\n'
            self._socket.send(cmd.encode('utf-8'))
            return self._socket.recv(4096).decode('utf-8').rstrip('\n')

    def save(self, path=None):
        if path is not None:
            # igv assumes the path is just a single filename, but
            # we can set the snapshot dir. then just use the filename.
            dirname = op.dirname(path)
            if dirname:
                self.set_path(dirname)
            return self.send('snapshot ' + op.basename(path))
        else:
            return self.send('snapshot')
    snapshot = save
def main_fn(igvsample,igvgene,bamfolder):

    print("BAM")
    #doctest.testmod()
    igv = IGV()
    igv.genome('hg38')
    bfname=igvsample
    print(bfname)
    #bamfolder="/mnt/seqdata/Exomsekvensering/Blod_exom/cnv/"
    bamfile=bamfolder+"{}-ready.bam".format(bfname)
    print(bamfile)
    igv.load(bamfile)
    #bfname=cnvtab4.sampleid
    print("genome")
    igv_gn=igvgene
    igv.go(igv_gn)
   
    print("bamloaded")
    
