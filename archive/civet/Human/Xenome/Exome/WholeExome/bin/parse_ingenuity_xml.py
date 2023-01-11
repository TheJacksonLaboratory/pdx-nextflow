#! /usr/bin/env python
"""
A quick script to parse the xml output by Ingenuity VCF, and turn it into tab-sep rows.
"""

import sys
import os
import argparse
from xml.dom import minidom

# cElementTree is faster, but not on all systems.
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

def parse_options():
    parser = argparse.ArgumentParser(description='LIMS Interface CGA pipelines')
    #parser.add_argument('--user', '-u')
    #parser.add_argument('--password', '-p')
    parser.add_argument('--output', '-o', help='filename for the output tabular text')
    parser.add_argument('--debug', action="store_true")
    # The input xml filename
    parser.add_argument('input_filename', help='input XML file name. Required.')
    args = parser.parse_args()
    return args

# DEBUGGING CODE    
DEBUG = True
def prettify_text(XML_text):
    reparsed = minidom.parseString(XML_text)
    return reparsed.toprettyxml(indent="  ")
    
def prettify(elem):
    rough_string = ET.tostring(elem, 'utf-8')
    return prettify_text(rough_string)

def main():
    args = parse_options()
    in_fn = args.input_filename
    if args.output:
        o_fn_base = args.output
    else:
        o_fn_base = in_fn
    o_fn_base = os.path.splitext(o_fn_base)[0]
    o_therapies = o_fn_base + '_therapies.txt'
    o_trials = o_fn_base + '_clinical_trials.txt'
    o_diagnoses = o_fn_base + '_diagnoses.txt'
    o_prognoses = o_fn_base + '_prognoses.txt'
    
    with open(in_fn) as f:
        tree = ET.parse(f)
    
    with open(o_therapies, 'w') as f:
        # Print the header
        print >> f, 'drug\ttype\tindication\treferences\tmatches\tgenes\tactivities\tpresent\ttranscript_id\ttranscript_hgvs\tprotein_id\tprotein_hgvs'
        for node in tree.findall('.//Therapy'):
            process_therapy(node, f)

    with open(o_trials, 'w') as f:
        print >> f, 'drug\tstudy\tphase\tstatus\tindication\treferences\tmatches\tgenes\tactivities\tpresent\ttranscript_id\ttranscript_hgvs\tprotein_id\tprotein_hgvs'
        for node in tree.findall('.//ClinicalTrial'):
            process_trial(node, f)

    with open(o_diagnoses, 'w') as f:
        print >> f, 'gene name'
        for node in tree.findall('.//Diagnosis'):
            process_diagnosis(node, f)

    with open(o_prognoses, 'w') as f:
        print >> f, 'gene_name'
        for node in tree.findall('.//Prognosis'):
            process_prognosis(node, f)

    print 'Done.'

def get_tag_text(node, tag):
    # Since we don't know whether any of the tags are listed multiple times,
    # assume they all are.
    l = []
    subs = node.findall('.//' + tag)
    for s in subs:
        l.append(s.text)
    return ', '.join(l)

def process_references(refs):
    l = []
    for ref in refs:
        src = ref.attrib['sourceId']
        typ = ref.attrib['sourceType']
        l.append('Source: {0}, Type: {1}'.format(src, typ))
    return ' | '.join(l)

def process_trial(trial, f):
    """
    <Drug>palbociclib</Drug>
    <Study>Phase 2 Trial of the Cyclin-Dependent Kinase Inhibitor PD 0332991 in Patients With Cancer</Study>
    <Phase>Phase 2</Phase>
    <Status>Recruiting</Status>
    <Indication>HER2-positive breast cancer</Indication>
    <references>
    <referenceInfo sourceId="NCT01037790" sourceType="NCT"></referenceInfo>
    </references>
    """
    drugs = get_tag_text(trial, 'Drug')
    studies = get_tag_text(trial, 'Study')
    phases = get_tag_text(trial, 'Phase')
    statuses = get_tag_text(trial, 'Status')
    indications = get_tag_text(trial, 'Indications')
    references = process_references(trial.findall('.//referenceInfo'))

    # Now process all the variants, emitting one line per variant.
    variants = trial.findall('.//VARIANT')
    for variant in variants:
        process_variant(variant, f, [drugs, studies, phases, statuses, indications, references])


def process_diagnosis(diagnosis, f):
    """
    
    At Guru's request, we're only going to output the gene names from this section of the input XML.
    
    NOTE: THE XML IN THIS SECTION IS !!!BADLY, HORRIBLY!!! BROKEN.  They have to fix it.  I'm not 
    going to do ANYTHING to try to clean it up.
    
    <Diagnosis>
    <references>
    <reference>
    <IngenuityFinding>&lt;pd&gt;&lt;show value=&quot;Mutant&amp;lt;b&amp;gt; human &quot;/&gt;&lt;show value=&quot;ASXL1&quot;&gt;&lt;hlink value=&quot;ING:2bcg&quot;/&gt;&lt;/show&gt;&lt;show value=&quot; gene is associated with myelodysplastic syndrome in&amp;lt;/b&amp;gt; human.&quot;/&gt;&lt;/pd&gt;</IngenuityFinding><referenceInfo sourceId="NCCN clinical Practice Guidelines in Oncology (NCCN Guidelines) for Myelodysplastic Syndromes V.2.2014" sourceType="URL"></referenceInfo>
    </reference>
    <reference>
    <IngenuityFinding>&lt;pd&gt;&lt;show value=&quot;Mutant&amp;lt;b&amp;gt; human &quot;/&gt;&lt;show value=&quot;ASXL1&quot;&gt;&lt;hlink value=&quot;ING:2bcg&quot;/&gt;&lt;/show&gt;&lt;show value=&quot; gene is associated with death of&amp;lt;/b&amp;gt;human exhibiting myelodysplastic syndrome.&quot;/&gt;&lt;/pd&gt;</IngenuityFinding><referenceInfo sourceId="NCCN clinical Practice Guidelines in Oncology (NCCN Guidelines) for Myelodysplastic Syndromes V.2.2014" sourceType="URL"></referenceInfo>

    Note: What they tried to do in the previous line is:
    <IngenuityFinding>
        <pd>
          <show value="Mutant<b> human "/>
          <show value="ASXL1">
            <hlink value="ING:2bcg"/>
          </show>
          <show value=" gene is associated with death of </b>human exhibiting myelodysplastic syndrome."/>
        </pd>
    </IngenuityFinding>
    <referenceInfo sourceId="NCCN clinical Practice Guidelines in Oncology (NCCN Guidelines) for Myelodysplastic Syndromes V.2.2014" sourceType="URL"></referenceInfo> 

    but note that they are burying <b> </b> into two separate literal values...  Bizarre.
    End of note...
    
    </reference></references>
    <VARIANT>
    <GENE_NAME>ASXL1</GENE_NAME>
    </VARIANT>
    </Diagnosis>
    """
    genes = []
    for variant in diagnosis.findall('.//VARIANT'):
        gene_names = get_tag_text(variant, 'GENE_NAME')
        for gene_name in gene_names.split(', '):
            if gene_name not in genes:
                genes.append(gene_name)
    for gene in genes:
        print >> f, gene

def process_prognosis(prognosis, f):
    """
    See header noted for process_diagnosis... Sigh.
    """
    genes = []
    for variant in prognosis.findall('.//VARIANT'):
        gene_names = get_tag_text(variant, 'GENE_NAME')
        for gene_name in gene_names.split(', '):
            if gene_name not in genes:
                genes.append(gene_name)
    for gene in genes:
        print >> f, gene

def process_therapy(therapy, f):
    """
    <Drug>trastuzumab</Drug>
    <Type>sensitive</Type>
    <Indication>inflammatory breast cancer</Indication>
    <references>
    <referenceInfo sourceId="NCCN clinical Practice Guidelines in Oncology (NCCN Guidelines) for Breast Cancer V.3.213" sourceType="URL"></referenceInfo>
    </references>
    """
    
    drugs = get_tag_text(therapy, 'Drug')
    types = get_tag_text(therapy, 'Type')
    indications = get_tag_text(therapy, 'Indication')
    references = process_references(therapy.findall('.//referenceInfo'))
    
    # Now process all the variants, emitting one line per variant.
    variants = therapy.findall('.//VARIANT')
    for variant in variants:
        process_variant(variant, f, [drugs, types, indications, references])

def process_variant(variant, f, others, ):
    """
    <VARIANT>
    <MATCH>gene</MATCH>
    <GENE_NAME>ERBB2</GENE_NAME>
    <GENE_ACTIVITY>gain</GENE_ACTIVITY>
    <PRESENT>true</PRESENT>
    <TRANSCRIPT>
    <TRANSCRIPT_ID>NM_004448.2</TRANSCRIPT_ID>
    <HGVS_TRANSCRIPT>c.1963A>G</HGVS_TRANSCRIPT>
    </TRANSCRIPT>
    <PROTEIN>
    <PROTEIN_ID>NP_004439.2</PROTEIN_ID>
    <HGVS_PROTEIN>p.I655V</HGVS_PROTEIN>
    </PROTEIN>
    </VARIANT>
    """
    matches = get_tag_text(variant, 'MATCH')
    genes = get_tag_text(variant, 'GENE_NAME')
    activities = get_tag_text(variant, 'GENE_ACTIVITY')
    presents = get_tag_text(variant, 'PRESENT')
    transcript_id, transcript_hgvs = process_transcripts(variant)
    protein_id, protein_hgvs = process_proteins(variant)

    # ... matches genes activities present transcript protein
    oths = '\t'.join(others)

    # Be careful!  The input XML is utf-8 format and can contain non-ASCII characters!
    print >> f, u'\t'.join([oths, matches, genes, activities, presents, 
        transcript_id, transcript_hgvs, protein_id, protein_hgvs]).encode('utf-8').strip()
        
def process_transcripts(node):
    transcripts = node.findall('.//TRANSCRIPT')
    # Guru believes that there is only one transcript per variant.
    # Emit an error message if there are more.
    if len(transcripts) == 0:
        return '', ''
    if len(transcripts) != 1:
        print >> sys.stderr, 'ERROR: There are more than one transcripts listed. Using the first.'
    xscr = transcripts[0]
    ids = get_tag_text(xscr, 'TRANSCRIPT_ID')
    hgvs = get_tag_text(xscr, 'HGVS_TRANSCRIPT')
    return ids, hgvs

def process_proteins(node):
    proteins = node.findall('.//PROTEIN')
    # Guru believes that there is only one protein per variant.
    # Emit an error message if there are more.
    if len(proteins) == 0:
        return '', ''
    if len(proteins) != 1:
        print >> sys.stderr, 'ERROR: There are more than one proteins listed. Using the first.'
    ptn = proteins[0]
    ids = get_tag_text(ptn, 'PROTEIN_ID')
    hgvs = get_tag_text(ptn, 'HGVS_PROTEIN')
    return ids, hgvs


    
main()