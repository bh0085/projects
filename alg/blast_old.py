import logging

from pylons import request, response, session, tmpl_context as c, url
from pylons.controllers.util import abort, redirect

from tal.lib.base import BaseController, render

log = logging.getLogger(__name__)

import tal.utils.security as sec
import subprocess as spc
import os, itertools as it
import Bio.Seq as bseq
import Bio.SeqIO as sio
import Bio.SeqRecord as sr
import Bio.Alphabet as alphabet

import json
from tal.model import User, Group, meta, Vector

import re

class BlastController(BaseController):

    def index(self):
        # Return a rendered template
        #return render('/blast.mako')
        # or, return a string
        return 'Hello World'

    def query(self):
        p = request.params;
        q = p.get('query',
                  'TATATATATGCGTTCTTTCATAGCTTGAATATCAAAAATGGGAAAATCTTGAAAAAAAATCCCAAAAAA'
                  )
        dogroups = p.get('groups', True)
        dousers = p.get('users', True)
        results = self._query(
            q,
            groups = dogroups,
            users = dousers)
        return results
    
    def _query(self,query,groups = True, users = True):
        
        all_users = meta.Session.query(User).all()
        all_groups = meta.Session.query(Group).all()

        results = {'groups':[],
                   'users':[]}

        hits = []
        for g in all_groups:
            root = '/data/blast/groups/{0}'.format(g.name)
            cmd = 'echo {0} | blast -p blastn -d db'\
                .format(query)
            prc = spc.Popen(cmd, 
                            stdout = spc.PIPE,
                            shell = True,
                            cwd = root)
            out = prc.stdout.read()
            
            reading = False;
            for l in out.splitlines():
                if len(l) == 0: continue;
                if l[0] == '>': 
                    reading = True;
                    entry = [];
                if not reading: continue;
                entry.append(l)
                if l[0:5] == 'Sbjct': 
                    reading = False
                    hits.append(entry)
    

        for u in all_users:
            root = '/data/blast/users/{0}'.format(u.username)
            cmd = 'echo {0} | blast -p blastn -d db'\
                .format(query)
            prc = spc.Popen(cmd, 
                            stdout = spc.PIPE,
                            shell = True,
                            cwd = root)
            out = prc.stdout.read()
            
            reading = False;
            for l in out.splitlines():
                if len(l) == 0: continue;
                if l[0] == '>': 
                    reading = True;
                    entry = [];
                if not reading: continue;
                entry.append(l)
                if l[0:5] == 'Sbjct': 
                    reading = False
                    hits.append(entry)
    
              

        hits_out = []
        for h in hits:
            l0 = h[0]
            match = re.compile('>([^_]+)_(\d+)').search(l0)
            owner_string = match.group(1)
            vector_id = match.group(2)
            matchlines = lines = re.compile('Query.*\n.*\n.*',re.M)\
                .search('\n'.join([l.strip() 
                                   for l in h])).group().splitlines()
        
            if 'group' in owner_string:
                owner_string = meta.Session.query(Group)\
                    .filter_by(id = owner_string[len('group'):])\
                    .first().name;
                
        
            if 'user' in owner_string:
                owner_string = 'User: '+meta.Session.query(User)\
                    .filter_by(uid = owner_string[len('user'):])\
                    .first().username;

            hits_out.append({
                    'text': '\n'.join(h),
                    'owner_string':owner_string,
                    'vector_id': vector_id,
                    'vector_name':meta.Session.query(Vector)\
                        .filter_by(id = vector_id).first().name,
                    'matchlines':matchlines
                    })
        return json.dumps(hits_out)
            

    
    def compile_all(self):
        for u in meta.Session.query(User).all():
            self.compile_userdb(u)
            
        for g in meta.Session.query(Group).all():
            self.compile_groupdb(g)

    
    def compile_userdb(self, user):
        root = '/data/blast/users/{0}'.format(user.username)
        if not os.path.isdir(root):
            os.makedirs(root)
        vecs = user.vectors;
        
        seqs, ids = [],[]
        for v in vecs:
            ids.append(v.id)
            seqs.append(v.record[0].sequence)
            
        fa_d = os.path.join(root, 'fasta')
        if not os.path.isdir(fa_d): os.mkdir(fa_d)
        

        fafiles = []
        for seq, idval in zip(seqs, ids):
            fname = os.path.join(fa_d,'user{1}_{0}.fa'\
                                     .format(idval,
                                             user.uid))
            fafiles.append(fname)
            sio.write(
                sr.SeqRecord(bseq.Seq(seq[0].seq, 
                                      alphabet = alphabet.DNAAlphabet), 
                             id = 'user{1}_{0}'.format(idval,
                                                       user.uid)
                             )
                ,fname,
                'fasta')
        cmd = 'formatdb -i "{0}" -p F -o T -n db'.format(' '.join(fafiles))
        prc = spc.Popen(cmd, 
                  shell = True,
                  cwd = root,
                  stdout = spc.PIPE)
        out =  prc.stdout.read()
        print out
        return out
        
            
    def compile_groupdb(self, group):
        root = '/data/blast/groups/{0}'.format(group.name)
        if not os.path.isdir(root):
            os.makedirs(root)

        users = [assoc.user for assoc in group.users]
        vecs =list(it.chain(*[user.vectors for user in users]))
        
        seqs, ids = [],[]
        for v in vecs:
            ids.append(v.id)
            seqs.append(v.record[0].sequence)
            
        fa_d = os.path.join(root, 'fasta')
        if not os.path.isdir(fa_d): os.mkdir(fa_d)
        

        fafiles = []
        for seq, idval in zip(seqs, ids):
            fname = os.path.join(fa_d,'group{1}_{0}.fa'\
                                     .format(idval,
                                             group.id))
            fafiles.append(fname)
            sio.write(
                sr.SeqRecord(bseq.Seq(seq[0].seq, 
                                      alphabet = alphabet.DNAAlphabet), 
                             id = 'group{1}_{0}'.format(idval,
                                                        group.id)
                             )
                ,fname,
                'fasta')
        cmd = 'formatdb -i "{0}" -p F -o T -n db'.format(' '.join(fafiles))
        prc = spc.Popen(cmd, 
                  shell = True,
                  cwd = root,
                  stdout = spc.PIPE)
        out =  prc.stdout.read()
        print out
        return out
        
        

    def search_user(self):
        user = sec.current_user();
        p = request.params;
        seq = p['sequence']
