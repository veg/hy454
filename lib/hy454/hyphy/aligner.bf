/*---------------------------------------------
 reverse complement a nucleotide string
---------------------------------------------*/

_revcomp_map = {};
_revcomp_map["A"] = "T";
_revcomp_map["C"] = "G";
_revcomp_map["G"] = "C";
_revcomp_map["T"] = "A";
_revcomp_map["M"] = "K";
_revcomp_map["R"] = "Y";
_revcomp_map["W"] = "W";
_revcomp_map["S"] = "S";
_revcomp_map["Y"] = "R";
_revcomp_map["K"] = "M";
_revcomp_map["B"] = "V";  /* not A */
_revcomp_map["D"] = "H";  /* not C */
_revcomp_map["H"] = "D";  /* not G */
_revcomp_map["V"] = "B";  /* not T */
_revcomp_map["N"] = "N";

function RevComp( _seq )
{
    _seqL = Abs( _seq );
    _seq2 = "";
    _seq2 * 128;
    _seqL += -1;
	for ( _rcidx = _seqL; _rcidx >= 0; _rcidx += -1 )
	{
		_seq2 * _revcomp_map[ _seq[ _rcidx ] && 1 ];
	}
	_seq2 * 0;
	return _seq2;
}

function Uppercase( _str )
{
    _upstr = _str && 1;
    _upstr * 0;
    return _upstr;
}

function TrimStop( _seq )
{
    _seq2 = _seq  ^ {{ "[-]$", "" }};
    _seq2 = _seq2 ^ {{ "[Tt][Aa][Aa]$", "" }};
    _seq2 = _seq2 ^ {{ "[Tt][Aa][Gg]$", "" }};
    _seq2 = _seq2 ^ {{ "[Tt][Gg][Aa]$", "" }};
    return _seq2;
}

function CleanAlignment( _aln, _keepIns )
/*
 * Given the raw alignment record from AlignSequence, clean the aligned sequence and return either the raw sequence or the list of positions by CSV
 * @param  _aln 		-- the alignment dictionary
 * @param  _keepIns 	-- whether or not to keep insertions in the REFERENCE (0/1)
 * @return the cleaned string
 */
{
    _ref = ( _aln[0] )[1];
    _seq = ( _aln[0] )[2];
    _overlap_count = 0;
    if ( _keepIns ) {
        _refLen = Abs( _ref );
        for ( _l = 0; _l < _refLen; _l += 1 ) {
            if ( _ref[ _l ] == _seq[ _l ] ) {
                _overlap_count += 1;
            }
        }
        _newRef = _ref;
        _newSeq = _seq;
    } else {
        _newRef = _ref ^ {{ "[-]",   ""  }};
        _newSeq = "";
        _newSeq * Abs( _newRef );

        _refLen = Abs( _ref );
        for ( _l = 0; _l < _refLen; _l += 1 ) {
            if ( _ref[ _l ] == _seq[ _l ] ) {
                _overlap_count += 1;
            }
            if ( _ref[ _l ] != "-" ) {
                _newSeq * (_seq[ _l ]); 
            }
        }
    }

    return { "ref": Uppercase( _newRef ),
             "seq": Uppercase( _newSeq ),
             "overlap": _overlap_count };
}

// -------------------------------------------------------------------------- //

function computeExpectedPerBaseScore( _expectedIdentity ) {
    meanScore = 0;

    for (_n1 = 0; _n1 < 4; _n1 += 1) {
        for (_n2 = 0; _n2 < 4; _n2 += 1) {
            if ( _n1 != _n2 ) {
                meanScore += ( 1 - _expectedIdentity ) * _aln_scorematrix[_n1][_n2] * _aln_base_freqs[_n1] * _aln_base_freqs[_n2] * _aln_pair_norm;
            } else {
                meanScore += _expectedIdentity * _aln_scorematrix[_n1][_n1] * _aln_base_freqs[_n1];
            }
        }
    }

    return meanScore;
}


// -------------------------------------------------------------------------- //
// ---------------------------- BEGIN MAIN ---------------------------------- //
// -------------------------------------------------------------------------- //



// _aln_scorematrix =
// {
//  { 5, -4, -4, -4}
//  {-4,  5, -4, -4}
//  {-4, -4,  5, -4}
//  {-4, -4, -4,  5}
// };

// 40% GC content
_aln_base_freqs = {
{0.30}
{0.20}
{0.20}
{0.30}
};

_aln_pair_norm = 0;
for ( _n1 = 0; _n1 < 4; _n1 += 1 ) {
    _aln_pair_norm += _aln_base_freqs[ _n1 ] * _aln_base_freqs[ _n1 ];
}
_aln_pair_norm = 1 / ( 1 - _aln_pair_norm );

ChoiceList ( _aln_dorevcomp, "Align reverse complement?", 1, SKIP_NONE,
             "No", "Do not check reverse complement",
             "Yes", "Check reverse complement" );

ChoiceList ( _aln_keepins, "Keep insertions?", 1, SKIP_NONE,
             "No", "Do not keep insertions",
             "Yes", "Keep insertions" );

// fprintf( stdout, "Please provide the reference sequence: " );
fscanf( stdin, "String", _aln_refseq );
// fprintf( stdout, "Please provide the expected identity: " );
fscanf( stdin, "Number", _aln_expected_identity );
// fprintf( stdout, "Please provide the number of query sequences: " );
fscanf( stdin, "Number", _aln_numseqs );
_aln_seqs = {};

for ( _aln_idx = 0; _aln_idx < _aln_numseqs; _aln_idx += 1 ) {
    // fprintf( stdout, "Please provide query sequence no. " + ( _aln_idx + 1 ) + ": " );
    fscanf( stdin, "String", _aln_grabseq );
    _aln_seqs[ _aln_idx ] = _aln_grabseq;
}

// uppercase and gapless
_aln_refseq = Uppercase( _aln_refseq ^ {{ "[-]", "" }} );
for ( _aln_idx = 0; _aln_idx < _aln_numseqs; _aln_idx += 1 ) {
    _aln_seqs[ _aln_idx ] = Uppercase( _aln_seqs[ _aln_idx ] ^ {{ "[-]", "" }} );
}

_aln_alnopts = {};
_aln_alnopts ["SEQ_ALIGN_SCORE_MATRIX"] = _aln_scorematrix;
_aln_alnopts ["SEQ_ALIGN_GAP_OPEN"] = 40;
_aln_alnopts ["SEQ_ALIGN_AFFINE"] = 1;
_aln_alnopts ["SEQ_ALIGN_GAP_OPEN2"] = 20;
_aln_alnopts ["SEQ_ALIGN_GAP_EXTEND2"] = 1;
_aln_alnopts ["SEQ_ALIGN_GAP_EXTEND"] = 10;
_aln_alnopts ["SEQ_ALIGN_NO_TP"] = 1; // this means local alignment, apparently
_aln_alnopts ["SEQ_ALIGN_CHARACTER_MAP"] = _aln_letters;

_aln_outstr = "";
_aln_outstr * 256;
_aln_outstr = "[";
// _aln_scores = { 1, Abs( _aln_seqs ) };

// if the expected_identity score is 0, then don't compute identity scores
_aln_expected_identity_score = 0;
if ( _aln_expected_identity > 0 ) {
    _aln_expected_identity_score = computeExpectedPerBaseScore( _aln_expected_identity );
}

// fprintf( stdout, "" + _aln_numseqs + "\n" );

for ( _aln_idx = 0; _aln_idx < _aln_numseqs; _aln_idx += 1 )
{
    // get the input sequence length, so we can normalize the score later
    _aln_seqlen = Abs( _aln_seqs[ _aln_idx ] );
    // align the sequences and get the score
    _aln_inseqs = {{ _aln_refseq, _aln_seqs[ _aln_idx ] }};
    AlignSequences ( _aln_alnseqs, _aln_inseqs, _aln_alnopts );
    // fprintf( stdout, "\n" + ( _aln_alnseqs[0] )[1] + "\n" + ( _aln_alnseqs[0] )[2] + "\n" );
    // divide the score by the length of the input sequence (which is assumed to be gapless)
    _aln_score = ( _aln_alnseqs[0] )[0] / _aln_seqlen;
    if ( _aln_dorevcomp ) {
        // if we are going to check the reverse complement:
        // once again align the sequences ( this time with the reverse complement )

        _aln_inseqs_rc = {{ _aln_refseq, RevComp( _aln_seqs[ _aln_idx ] ) }};
        AlignSequences( _aln_alnseqs_rc, _aln_inseqs_rc, _aln_alnopts );
        // divide the score by the length of the input sequence (which is assumed to be gapless)
        _aln_score_rc = ( _aln_alnseqs_rc[0] )[0] / _aln_seqlen;
        if ( _aln_score_rc > _aln_score ) {
            // if the reverse complement score is greater than the regular score, use it instead
            _aln_cleanseqs = CleanAlignment( _aln_alnseqs_rc, _aln_keepins );
            _aln_score = _aln_score_rc;
        } else {
            // otherwise just the regular score
            _aln_cleanseqs = CleanAlignment( _aln_alnseqs, _aln_keepins );
        }
    } else {
        // if we are not checking the reverse complement, just score the result
        _aln_cleanseqs = CleanAlignment( _aln_alnseqs, _aln_keepins );
    }

    // trim the sequence back to the true length (unpadded) if we're not keeping insertions
    // and modify the alignment score to account for this
    _aln_cleanseq = _aln_cleanseqs[ "seq" ];

    if (_aln_idx > 0)
    {
        _aln_outstr * ",";
    }
    _aln_identity_score = _aln_score - _aln_expected_identity_score;
    _aln_outstr * ( "[\"" + _aln_cleanseq + "\"," + _aln_score + "," + _aln_cleanseqs["overlap"] + "," + _aln_identity_score + "]" );
}

_aln_outstr * "]";
_aln_outstr * 0;

