#!/usr/bin/python

import sys

frm_usage = '<p style="margin-left: 0.19in; margin-bottom: 0.08in"><font face="monospace" size="2"><b>{}</b></font></p>'

frm_option = '<tr valign="top"><td width="34%"><p style="margin-left: 0.39in; margin-bottom: 0.0cm; margin-top: 0.0cm">' \
  '<font face="monospace" size="2"><b>{}</b></font></p></td>' \
  '<td width="66%"><p style="margin-left: 0.39in; margin-bottom: 0.0cm; margin-top: 0.0cm">' \
  '<font face="monospace" size="2"><b>{}</b></font></p></td></tr>'

hdr = '<table width="100%" border="0" cellpadding="0" cellspacing="0"><col width="88*"><col width="168*"><tbody>'
ftr = '</tbody></table>'


if __name__ == "__main__":
    intab = False
    for l in sys.stdin:

        if l.startswith("Usage: "):
            if intab:
                print(ftr)
                intab = False
            le = l[7:].replace('<', '&lt;').replace('>', '&gt;')
            print(frm_usage.format(le))
        else:
            if not intab:
                print(hdr)
                intab = True
            kw = l[:23].strip(' ').replace('<', '&lt;').replace('>', '&gt;').replace(' ', '&nbsp;')
            desc = l[23:].replace('<', '&lt;').replace('>', '&gt;')
            print(frm_option.format(kw, desc))

    if intab:
        print(ftr)
