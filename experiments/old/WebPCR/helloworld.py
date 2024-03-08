#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cgi
import webapp2
from google.appengine.api import users

import string
from pcr import Anneal
from parse import parse


template_separator  =u'templates'
flankuplength       = 200
flankdnlength       = 200
max_product_size    = 10000
homology_limit      = 12
cutoff_featured_template = 6
cutoff_detailed_figure   = 5
cutoff_minimum_lines_for_report_in_subpage = 20

report_header=u'''
==============
PCR simulation
==============

'''

report_for_each_simulation = u'$anneal_primers'

report_for_each_amplicon =  u'''
PCR product from $template_name:

>$product_name
$product_sequence

$figure
suggested PCR programs
$pcr_program
'''

#uncomment the following lines to get flanking sequences in report
#>upstream_flanking_sequence
#$upstream_flanking_sequence
#>downstream_flanking_sequence
#$downstream_flanking_sequence
#---

#settings for flavios primers
#report_header=""
#report_for_each_simulation = ""
#report_for_each_amplicon = '''$forward_primer_name
#$reverse_primer_name
#$figure
#>upstream_flanking_sequence
#$upstream_flanking_sequence
#---


class MainPage(webapp2.RequestHandler):
    def get(self):
        self.response.out.write("""
        WebPCR
          <html>
           <head>
            <title>WebPCR</title> version 0.50<br>
            </head>

            <br> This service simulate PCR given a list of sequences in FASTA or Genbank format.
            <br> The last sequence in the list is interpreted as template and the preceding ones as primers.
            <br> For each potential PCR product, a report will be generated containing the PCR product sequence and some additional data.
            <br> Delete the example content below and past your own data. <a href="https://docs.google.com/document/d/1JgmwWSL6axFw3Ed9aKurpC_qwo_PoILAxvpbUU5nviY/edit?usp=sharing">More help</a>.
            <br><br>


            <body>
              <form action="/sign" method="post">
                <div><textarea name="content" rows="20" cols="120">
>ForwardPrimer
gctactacacacgtactgactg

>ReversePrimer
tgtggttactgactctatcttg

>MyTemplate
gctactacacacgtactgactgcctccaagatagagtcagtaaccaca
                </textarea></div>
                <div><input type="submit" value="Run simulation"></div>
              </form>
            </body>
          </html>""")
class Guestbook(webapp2.RequestHandler):
    def post(self):
        self.response.out.write('<html><body><pre>')
        input = self.request.get('content')


        seqs = parse(input)

        if not seqs:
            return

        template = seqs.pop()



        result_text=''
        message_template = report_header

        anneal_primers = Anneal( seqs,
                                 template,
                                 homology_limit,
                                 max_product_size)

        if anneal_primers.number_of_products==0:
            result_text="\n"+anneal_primers.report()

        elif 1<=anneal_primers.number_of_products<=cutoff_detailed_figure:
            message_template += report_for_each_simulation
            for amplicon in anneal_primers.amplicons:
                message_template += report_for_each_amplicon
                result_text+="\n"+string.Template(message_template).safe_substitute(
                    anneal_primers                = anneal_primers,
                    forward_primer_name           = amplicon.forward_primer.primer.name,
                    forward_primer_sequence       = amplicon.forward_primer.primer.seq,
                    reverse_primer_name           = amplicon.reverse_primer.primer.name,
                    reverse_primer_sequence       = amplicon.reverse_primer.primer.seq,
                    product_name                  = amplicon.pcr_product().id,
                    product_sequence              = amplicon.pcr_product().seq,
                    template_name                 = anneal_primers.template.name,
                    template_sequence             = anneal_primers.template.seq,
                    upstream_flanking_sequence    = amplicon.flankup(),
                    downstream_flanking_sequence  = amplicon.flankdn(),
                    figure                        = amplicon.detailed_figure(),
                    pcr_program                   = amplicon.pcr_program() )
                message_template=''

        elif anneal_primers.number_of_products>cutoff_featured_template:
            result_text+="\n"+anneal_primers.featured_template().format("gb")


        self.response.out.write(cgi.escape(result_text))
        self.response.out.write('</pre></body></html>')

app = webapp2.WSGIApplication([('/', MainPage),
                              ('/sign', Guestbook)],
                              debug=True)
