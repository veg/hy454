
from os.path import abspath, exists, join, split
from sys import stderr
from itertools import product
from operator import itemgetter
import json

__all__ = ['DRAM','codon_mapper_max', 'codon_mapper_min', 'codon_mapper_avg']

_codonTable = {'AAA' : 'K', 'AAC' : 'N', 'AAG' : 'K', 'AAT' : 'N', 'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T', 'AGA' : 'R', 'AGC' : 'S', 'AGG' : 'R', 'AGT' : 'S', 'ATA' : 'I', 'ATC' : 'I', 'ATG' : 'M', 'ATT' : 'I', 'CAA' : 'Q', 'CAC' : 'H', 'CAG' : 'Q', 'CAT' : 'H', 'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P', 'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R', 'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L', 'GAA' : 'E', 'GAC' : 'D', 'GAG' : 'E', 'GAT' : 'D', 'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A', 'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G', 'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V', 'TAA' : 'X', 'TAC' : 'Y', 'TAG' : 'X', 'TAT' : 'Y', 'TCA' : 'S', 'TCC' : 'S', 'TCG' : 'S', 'TCT' : 'S', 'TGA' : 'X', 'TGC' : 'C', 'TGG' : 'W', 'TGT' : 'C', 'TTA' : 'L', 'TTC' : 'F', 'TTG' : 'L', 'TTT' : 'F'}
_nucCodes   = {'A' : ['A'], 'C' : ['C'], 'G' : ['G'], 'T' : ['T'], 'M' : ['A','C'], 'R' : ['A', 'G'], 'W' : ['A', 'T'], 'S' : ['C', 'G'], 'Y' : ['C','T'], 'K' : ['G', 'T'], 'V': ['A','C','G'], 'H' : ['A','C','T'], 'D' : ['A','G','T'], 'B' : ['C', 'G', 'T'], 'N': ['A','C','G','T']}

_resistance_files = {'PI' : {'name': 'PI.dat', 'offset': 0, 'gene': 'PR'}, 'NRTI' : {'name': 'NRTI.dat', 'offset' : 99, 'gene': 'RT'}, 'NNRTI' : {'name': 'NNRTI.dat', 'offset' : 99, 'gene': 'RT'}}
_resistance_data  = None

COMPLETE_LIST, BINARY  = range(2)

def read_single_file (file, drug_class, offset, gene, positions):
    resistance_info = {}
    header_line =  file.readline().rstrip().split ('\t')

    assert len (header_line) > 2, "Expected at least 3 columns in the drug resistance file header line"
    assert header_line[0].upper() == "POSITION", "The first column is expected to be named Position"
    assert header_line[1].upper() == "AA", "The second column is expected to be named AA"

    drug_names = header_line[2:]
    positions_by_class = set ()
    positions_by_drug  = {}

    for drug in drug_names:
        positions_by_drug [drug] = set()

    for line in file.readlines():
        line_tabs = line.rstrip().split('\t')
        line_item = {}
        coord = int(line_tabs[0]) + offset -1
        for idx,drug in enumerate(drug_names):
            if int (line_tabs[2+idx]) != 0:
                line_item [drug] = int (line_tabs[2+idx])
                positions_by_drug [drug].add (coord)

        if len (line_item) == 0:
            continue

        line_item ['Max'] = max ([line_item[k] for k in drug_names if k in line_item])
        positions.add (coord)
        positions_by_class.add(coord)
        resistance_info [str(coord) + line_tabs[1]] =  line_item
        line_item ['Original'] = gene + line_tabs[0] + line_tabs[1]

    for drug in drug_names:
        positions_by_drug [drug] = list(positions_by_drug [drug])

    return (resistance_info,{'names':drug_names, 'positions':list(positions_by_class)},positions_by_drug)

def import_mdr_scores_from_stanford ():
    resistance_data_here = {}
    unique_positions     = set ()

    by_class = {}
    by_drug  = {}

    for drug_class in _resistance_files:
        with open (join(split(abspath(__file__))[0], 'data', 'scores', _resistance_files[drug_class]['name'])) as input_file:
            class_info = read_single_file (input_file, drug_class, _resistance_files[drug_class]['offset'], _resistance_files[drug_class]['gene'], unique_positions)
            resistance_data_here.update(class_info[0])
            by_class[drug_class] = class_info[1]
            by_drug.update(class_info[2])

    return {'mutations': resistance_data_here, 'positions': list(unique_positions), 'by class' : by_class, 'by drug': by_drug}



def translate_codon(codon):
    if codon in _codonTable:
        return [_codonTable[codon]]
    else:
        try:
            position_translations = [_nucCodes[letter] for letter in codon]
            amino_acids = set()
            for resolution in product (position_translations[0], position_translations[1], position_translations[2]):
                amino_acids.add(_codonTable[''.join(resolution)])

            return list(amino_acids)

        except KeyError:
            return []

def codon_mapper_raw (position,codon,resistance_table):
    codon_mapped = translate_codon (codon)
    line_item   = set()
    if len (codon_mapped) > 4: #too many degeneracies
        line_item.add ('Noisy')
    for codon_resolution in codon_mapped:
        try_code = str(position) + codon_resolution
        if try_code in resistance_table['mutations']:
        	line_item.add (try_code)
        else:
            line_item.add ('WT')
    return list(line_item)


def codon_mapper_max (raw_info):
    pass


def DRAM_calculator (positional_data, resistance_table = None):
# positional_data should be of the form
# {'seqid1':[[coord1,codon1],...] } where coord a 0-based CODON (e.g. in nucleotide space)
# coordinate and codon is the three letter codon

    if resistance_table is None:
        resistance_table = import_mdr_scores_from_stanford ()

    resistance_report = {}
    all_positions = set (resistance_table['positions'])

    for seq_id in positional_data:
        position_list = positional_data [seq_id]
        sequenced_positions = set ([ k[0] // 3 for k in position_list ]).intersection (all_positions)
        sequenced_codons    = {}

        this_sequence_report = {}

        for codon in position_list:
            if codon[0] // 3 in sequenced_positions:
                sequenced_codons [ codon[0] // 3 ] = codon[1]

        for a_codon in sequenced_codons:
            codon_map = codon_mapper_raw (a_codon, sequenced_codons[a_codon], resistance_table)
            if len (codon_map):
                this_sequence_report [str(a_codon)] = codon_map
            else:
                this_sequence_report [str(a_codon)] = ["WT"]

        not_sequenced = all_positions.difference (sequenced_positions)
        for a_codon in not_sequenced:
            this_sequence_report [str(a_codon)] = None

        resistance_report [seq_id] = this_sequence_report

    resistance_report ['>resistance table'] = resistance_table
    return resistance_report


def _DRAM_resolve_score (scores, ambiguity_resolution):
    if ambiguity_resolution == DRAM_AMBIGS_OPTIMISTIC:
        return min (scores.items(), key = itemgetter (1))
    else:
        return max (scores.items(), key = itemgetter (1))



DRAM_AMBIGS_PESSIMISTIC, DRAM_AMBIGS_OPTIMISTIC = range (2)
DRAM_REPORT_COMPLETE, DRAM_REPORT_BY_DRUG, DRAM_REPORT_BY_MUTATION, DRAM_REPORT_BY_CLASS = range (4)

def DRAM_reporter (positional_resistance_data, resistance_table = None, ambiguity_resolution = None, report_mode = None, mut_score_limit = 30, drug_resistance_limit = 30):
# positional_resistance_data is returned by DRAM_calculator
    if resistance_table is None:
        if '>resistance table' in positional_resistance_data:
            resistance_table = positional_resistance_data['>resistance table']
        else:
            resistance_table = import_mdr_scores_from_stanford ()
    if ambiguity_resolution is None:
        ambiguity_resolution = DRAM_AMBIGS_PESSIMISTIC
    if report_mode is None:
        report_mode = DRAM_REPORT_COMPLETE

    dram_report = {}

    # generate the full report by default, then condense if necessary
    # score by drug (with comments if sequenced completely or has noisy data): min and max
    # resistance mutations by drug class (can be limited to major)
    # resistance to classes e.g. NRTI : ['susceptible', 'resistant']


    for sequence_id in positional_resistance_data:
        if sequence_id == '>resistance table':
            continue
        sequence_info = positional_resistance_data[sequence_id]
        sequence_report = {'by drug' : {}}
        for drug_name in resistance_table['by drug']:
            total_score = 0
            covered_positions = 0
            noisy_positions   = 0
            missing_positions = 0
            major_DRAM        = []

            list_of_positions_by_drug = resistance_table['by drug'][drug_name]

            for position in list_of_positions_by_drug:
                score = 0
                in_this_sequence = sequence_info[str(position)]
                resolution_scores = {}

                if in_this_sequence is None:
                    missing_positions += 1
                    continue
                elif in_this_sequence[0] == 'Noisy':
                    noisy_positions += 1
                else:
                    covered_positions += 1
                    for value in in_this_sequence:
                        if value != 'WT':
                            try:
                                score = resistance_table['mutations'][value][drug_name]
                            except KeyError:
                                score = 0
                        else:
                            score = 0
                        resolution_scores [value] = score
                    mutation,score  = _DRAM_resolve_score (resolution_scores, ambiguity_resolution)
                    if score >= mut_score_limit:
                        major_DRAM.append(resistance_table["mutations"][mutation]["Original"])

                total_score += score

            sequence_report ['by drug'][drug_name] = {'covered' : covered_positions, 'missing' : missing_positions, 'noisy': noisy_positions, 'score':  total_score ,'mutations' : major_DRAM if len (major_DRAM) else None}


        sequence_report ['by class'] = {}
        for drug_class, class_info in resistance_table['by class'].items():
            report = {'mutations' : set(), 'resistant' : []}
            for drug_id in class_info['names']:
                #print (drug_id, report)
                if sequence_report ['by drug'][drug_id]['score'] >= drug_resistance_limit:
                    report['resistant'].append (drug_id)
                if sequence_report ['by drug'][drug_id]['mutations'] is not None:
                    report['mutations'].update (sequence_report ['by drug'][drug_id]['mutations'])
            report['mutations'] = list (report['mutations'])
            sequence_report ['by class'][drug_class] = report

        dram_report [sequence_id] = sequence_report
        #assert False

    if report_mode == DRAM_REPORT_BY_DRUG:
    # in this mode the report is of the form
    # seq id => list of drugs to which the sequence is resistant, grouped by class
        drug_report = {}
        for seq_id in dram_report:
            seq_drug_report = {}
            for drug_class, drug_value in dram_report[seq_id]['by class'].items():
                seq_drug_report [drug_class] = drug_value['resistant']
            drug_report [seq_id] = seq_drug_report
        return drug_report
    elif report_mode == DRAM_REPORT_BY_MUTATION:
    # in this mode the report is of the form
    # seq id => list of major mutations
        mut_report = {}
        for seq_id in dram_report:
            seq_mut_report = set()
            for drug_class, drug_value in dram_report[seq_id]['by class'].items():
                seq_mut_report.update(drug_value['mutations'])
            mut_report [seq_id] = list(seq_mut_report)
        return mut_report
    elif report_mode == DRAM_REPORT_BY_CLASS:
    # in this mode the report is of the form
    # seq id => list of drug classes to which there are major mutations
        class_report = {}
        for seq_id in dram_report:
            class_list = set()
            for drug_class, drug_value in dram_report[seq_id]['by class'].items():
                if len (drug_value['resistant']):
                    class_list.add (drug_class)
            class_report [seq_id] = list(class_list)
        return class_report


    return dram_report

if __name__ == "__main__":
    from sys import argv
    assert len (argv) > 1, "This program requires a single argument -- the .json file of sequences split by positions"

    mode = None

    if len (argv) == 2:
        with open (argv[1]) as positional_output:
            DRAM_data = DRAM_calculator (json.load (positional_output))
    else:
        with open (argv[1]) as positional_output:
            DRAM_data = json.load (positional_output)
            mode = int (argv[2])

    print (json.dumps(DRAM_reporter (DRAM_data,report_mode = mode),sort_keys=True, indent=4))

