#!/usr/bin/env python
import os
import subprocess as spc

colors = ['white','lightgray', 'black', 'darkgray']
sizes = [128, 64, 32]

def main(make_imgs = False):
    print [os.path.splitext(s) for s in os.listdir(os.getcwd())]
    svgs = [ s for s in os.listdir(os.getcwd()) if os.path.splitext(s)[-1] == '.svg']
    svgs = [ s for s in svgs if not 'color=' in s]
    print svgs
    for s in svgs:
        for c in colors:
            path = os.path.join(os.getcwd(), s)
            prc = spc.Popen( 'svg_cvcolors -c {0} -a'.format(c),
                             stdin = spc.PIPE,
                             stdout = spc.PIPE,
                             shell = True)


            color_file =   '.'.join(os.path.splitext(s)[:-1]) \
                                             + '_color={0}'.format(c)
            color_outpath = os.path.join(os.getcwd(), color_file) + '.svg'
            
            inp_open = open(path)
            data = inp_open.read()
            out = prc.communicate(input = data)[0]
            inp_open.close()

            out_open = open(color_outpath, 'w')
            out_open.write(out)
            out_open.close()
            
            if make_imgs:
                for r in sizes:
                    size_file = color_file + '_size={0}'.format(r)
                    size_outpath = os.path.join(os.getcwd(),size_file + '.png')
                                            
                    cmd = 'rsvg -h {0} -w {0} {1} {2}'\
                        .format(r, color_outpath, size_outpath)
                    spc.Popen(cmd,shell = True).communicate()
            

    
if __name__ == "__main__":
    imgs = False
    main(make_imgs = imgs)
