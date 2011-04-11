#!/usr/bin/env python
'''
Install anNote.
'''
import os
import mailman

att_dir = os.path.join(os.environ.get('HOME'),'.anNote')
if not os.path.isdir(att_dir):
  print 'installing anNote...'
  print '  making ~/.anNote'
  os.mkdir(att_dir)
