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
    _newRef = _ref ^ {{ "[-]",   ""  }};
    _altRef = _ref ^ {{ "[a-z]", "_" }};
    _newStr = "";
    if ( _keepIns ) {
        _keepIns = 1;
        _newStr * Abs( _seq );
    } else {
        _keepIns = 0;
        _newStr * Abs( _newRef );
    }
    _k = 0;
    _overlap_count = 0;
    
    // codon by codon...
    _newRefLen = Abs( _newRef );
    for ( _l = 0; _l < _newRefLen; _l += 3 ) {
        _ins = 0;
        _dels = 0;
        // count the number of insertions and deletions
        _cdnAdv = Min( 3, _newRefLen - _l );
        for ( _k2 = 0; _k2 < _cdnAdv; _k2 += 1 ) {
            if (_altRef[ _k+_k2 ] == "-") {
                _ins += 1;
            } else {
                if ( _altRef[ _k+_k2 ] == "_" ) {
                    _dels += 1;
                }
            }
        }
        // if _ins == 3 (and we want to keep inserts),
        // then we have a full codon insertion,
        // add the original characters back in
        // if _ins == 0 and _dels == 0, then everything is normal,
        // add the original characters back in
        if ( ( _keepIns * _ins ) == 3 || ( _ins == 0 && _dels == 0 ) ) {
            _newStr * _seq[ _k ][ _k+2 ];
            _k += 3;
            
            _overlap_count += 3*( _seq[ _k ][ _k+2 ] != "---" );
        }
        // if neither of those two cases is true, then we need to go
        // position by position, removing insertions and
        // fixing deletions by adding in a "N"
        else {
            // _k2 advances by only a single codon (3 positions), whereas
            // _l2 moves us ahead in the alignment. At the end, the difference
            // between _k2 and _l2 should be that _l2 is greater by the number
            // of inserted nucleotides
            _k2 = 0;
            for ( _l2 = 0; _k2 < _cdnAdv; _l2 += 1 ) {
                // "_" means a deletion, add an "N" back in
                if ( _altRef[ _k+_l2 ] == "_" ) {
                    _newStr * "N";
                    _k2 += 1;
                }
                // if the character isn't an insertion
                // which it always isn't because we took care of
                // insertions above, add the original back in
                else {
                    if ( _altRef[ _k+_l2 ] != "-" ) {
                        _newStr * _seq[ _k+_l2 ];
                        _k2 += 1;
                        _overlap_count += (_seq[ _k+_l2 ] != "-");
                    }
                }
            }
            _k += _l2;
        }
    }
    _newRef * 0;
    _newStr * 0;
    // get rid of any gaps
    // _newStr2 = _newStr^{{"[-]", ""}};
    return { "ref": Uppercase( _newRef ),
             "seq": Uppercase( _newStr[0][ _newRefLen ] ),
             "overlap": _overlap_count/3};
}

function pSM2cSM(_protScoreMatrix, _protLetters)
{
    LoadFunctionLibrary ("chooseGeneticCode", {"00":"Universal"});
    LoadFunctionLibrary ("GrabBag");

    _scoreMatrix  = {65,65};
    _mapping      = mapStrings (_hyphyAAOrdering, _protLetters);

    for (_k = 0; _k < 64; _k += 1)
    {
        _mappedK = _mapping[_Genetic_Code[_k]];
        if (_mappedK >= 0)
        {
            for (_k2 = _k; _k2 < 64; _k2 += 1)
            {
                _mappedK2 = _mapping[_Genetic_Code[_k2]];
                if (_mappedK2 >= 0)
                {
                    _aScore = _protScoreMatrix[_mappedK][_mappedK2];
                    if (_mappedK == _mappedK2 && _k2 > _k)
                    {
                        _aScore = _aScore - 1;
                    }
                }
                else
                {
                    _aScore = -10000;
                }
                _scoreMatrix[_k][_k2] = _aScore;
                _scoreMatrix[_k2][_k] = _aScore;
            }
        }
        else
        {
            for (_k2 = _k; _k2 < 64; _k2 += 1)
            {
                _scoreMatrix[_k][_k2] = -10000;
                _scoreMatrix[_k2][_k] = -10000;
            }
        }
    }
    return _scoreMatrix;
}

function cSM2partialSMs(_scoreMatrix)
{
    m3x2  =  {65,48};
    m3x1  =  {65,12};

    for (thisCodon = 0; thisCodon < 64; thisCodon += 1)
    {
        for (d1 = 0; d1 < 4; d1 += 1)
        {
            max100 = -1e100;
            max010 = -1e100;
            max001 = -1e100;

            for (d2 = 0; d2 < 4; d2 += 1)
            {
                partialCodon = 4*d1 + d2;
                t = 16*d1 + d2;
                max110 = -1e100;
                max101 = -1e100;
                max011 = -1e100;

                for (d3 = 0; d3 < 4; d3 += 1)
                {
                    // d1 is 1
                    max100 =  Max(max100,_scoreMatrix[thisCodon][d1*16+d2*4+d3]);
                    max010 =  Max(max010,_scoreMatrix[thisCodon][d2*16+d1*4+d3]);
                    max001 =  Max(max001,_scoreMatrix[thisCodon][d2*16+d3*4+d1]);

                    // d1 and d2 are 1
                    max110 = Max(max110,_scoreMatrix[thisCodon][4*partialCodon + d3]);
                    max101 = Max(max101,_scoreMatrix[thisCodon][t + 4*d3]);
                    max011 = Max(max011,_scoreMatrix[thisCodon][partialCodon + 16*d3]);
                }
                m3x2[thisCodon][3*partialCodon]   = max110;
                m3x2[thisCodon][3*partialCodon+1] = max101;
                m3x2[thisCodon][3*partialCodon+2] = max011;
            }
            m3x1[thisCodon][3*d1]   = max100;
            m3x1[thisCodon][3*d1+1] = max010;
            m3x1[thisCodon][3*d1+2] = max001;
        }
    }
    return {"3x1": m3x1, "3x2": m3x2};
}

// -------------------------------------------------------------------------- //

function computeExcpectedPerBaseScore () {
    meanScore = 0;
    
    for (_aa1 = 0; _aa1 < 20; _aa1 += 1) {
        for (_aa2 = 0; _aa2 < 20; _aa2 += 1) {
            meanScore += _cdnaln_protScoreMatrix[_aa1][_aa2] * _cdaln_base_frequencies[_aa1] * _cdaln_base_frequencies[_aa2];
        }
    }
    
    return meanScore;
}


// -------------------------------------------------------------------------- //
// ---------------------------- BEGIN MAIN ---------------------------------- //
// -------------------------------------------------------------------------- //



_cdnaln_protScoreMatrix =
{
 { 6, -3, -4, -4, -2, -2, -2, -1, -3, -3, -3, -2, -2, -4, -2,  0, -1, -5, -3, -1, -4, -2, -2, -7}
 {-3,  8, -2, -4, -6,  0, -2, -5, -2, -6, -4,  1, -3, -5, -4, -2, -3, -5, -3, -5, -3, -1, -2, -7}
 {-4, -2,  8,  0, -5, -1, -2, -2,  0, -6, -6, -1, -4, -5, -4,  0, -1, -7, -4, -5,  6, -2, -2, -7}
 {-4, -4,  0,  8, -6, -2,  0, -3, -3, -5, -6, -2, -6, -6, -3, -1, -3, -7, -6, -6,  6,  0, -3, -7}
 {-2, -6, -5, -6, 10, -5, -7, -5, -5, -3, -3, -6, -3, -5, -5, -2, -2, -4, -4, -2, -5, -6, -4, -7}
 {-2,  0, -1, -2, -5,  8,  1, -4,  0, -6, -4,  0, -1, -6, -3, -1, -2, -3, -3, -4, -1,  6, -2, -7}
 {-2, -2, -2,  0, -7,  1,  7, -4, -1, -6, -5,  0, -4, -6, -3, -1, -2, -5, -4, -4,  0,  6, -2, -7}
 {-1, -5, -2, -3, -5, -4, -4,  7, -4, -7, -6, -3, -5, -5, -4, -2, -4, -4, -5, -6, -2, -4, -4, -7}
 {-3, -2,  0, -3, -5,  0, -1, -4, 10, -6, -5, -2, -3, -3, -4, -2, -4, -5,  0, -6, -1, -1, -3, -7}
 {-3, -6, -6, -5, -3, -6, -6, -7, -6,  6,  0, -5,  0, -1, -5, -5, -2, -5, -3,  2, -5, -6, -2, -7}
 {-3, -4, -6, -6, -3, -4, -5, -6, -5,  0,  6, -5,  1, -1, -5, -5, -3, -3, -3,  0, -6, -5, -2, -7}
 {-2,  1, -1, -2, -6,  0,  0, -3, -2, -5, -5,  7, -3, -6, -2, -1, -2, -5, -3, -4, -2,  0, -2, -7}
 {-2, -3, -4, -6, -3, -1, -4, -5, -3,  0,  1, -3,  9, -1, -5, -3, -2, -3, -3,  0, -5, -2, -1, -7}
 {-4, -5, -5, -6, -5, -6, -6, -5, -3, -1, -1, -6, -1,  8, -6, -4, -4,  0,  1, -3, -6, -6, -3, -7}
 {-2, -4, -4, -3, -5, -3, -3, -4, -4, -5, -5, -2, -5, -6,  9, -2, -3, -6, -5, -4, -4, -3, -4, -7}
 { 0, -2,  0, -1, -2, -1, -1, -2, -2, -5, -5, -1, -3, -4, -2,  7,  0, -5, -3, -4, -1, -1, -2, -7}
 {-1, -3, -1, -3, -2, -2, -2, -4, -4, -2, -3, -2, -2, -4, -3,  0,  7, -4, -3, -1, -2, -2, -2, -7}
 {-5, -5, -7, -7, -4, -3, -5, -4, -5, -5, -3, -5, -3,  0, -6, -5, -4, 12,  0, -6, -7, -4, -4, -7}
 {-3, -3, -4, -6, -4, -3, -4, -5,  0, -3, -3, -3, -3,  1, -5, -3, -3,  0,  9, -3, -5, -3, -2, -7}
 {-1, -5, -5, -6, -2, -4, -4, -6, -6,  2,  0, -4,  0, -3, -4, -4, -1, -6, -3,  6, -6, -4, -2, -7}
 {-4, -3,  6,  6, -5, -1,  0, -2, -1, -5, -6, -2, -5, -6, -4, -1, -2, -7, -5, -6,  7, -1, -3, -7}
 {-2, -1, -2,  0, -6,  6,  6, -4, -1, -6, -5,  0, -2, -6, -3, -1, -2, -4, -3, -4, -1,  7, -2, -7}
 {-2, -2, -2, -3, -4, -2, -2, -4, -3, -2, -2, -2, -1, -3, -4, -2, -2, -4, -2, -2, -3, -2, -2, -7}
 {-7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7,  1}
};

_cdaln_base_frequencies = {
{   0.060490222}
{   0.020075899}
{   0.042109048}
{   0.071567447}
{   0.028809447}
{   0.072308239}
{   0.022293943}
{   0.069730629}
{   0.056968211}
{   0.098851122}
{   0.019768318}
{   0.044127815}
{   0.046025282}
{   0.053606488}
{   0.066039665}
{    0.05060433}
{   0.053636813}
{   0.061625237}
{   0.033011601}
{   0.028350243}
};


ChoiceList ( _cdnaln_dorevcomp, "Align reverse complement?", 1, SKIP_NONE,
             "No", "Do not check reverse complement",
             "Yes", "Check reverse complement" );
             
fscanf( stdin, "String", _cdnaln_refseq );
fscanf( stdin, "Number", _cdnaln_numseqs );
_cdnaln_seqs = {};

for ( _cdnaln_idx = 0; _cdnaln_idx < _cdnaln_numseqs; _cdnaln_idx += 1 ) {
    fscanf( stdin, "String", _cdnaln_grabseq );
    _cdnaln_seqs[ _cdnaln_idx ] = _cdnaln_grabseq;
}

// Due to some bugs in the implementation of the codon aligner,
// pad the reference sequence to a multiple of 3, 
_cdnaln_truelen = Abs( _cdnaln_refseq );
_cdnaln_pad = 3 - ( _cdnaln_truelen % 3 );
if ( _cdnaln_pad > 0 && _cdnaln_pad < 3 ) {
    _cdnaln_refseq * ( "NN"[ 0 ][ ( _cdnaln_pad - 1 ) ] );
}
_cdnaln_refseq * 0;
_cdnaln_scoremod = _cdnaln_truelen / ( _cdnaln_truelen + _cdnaln_pad );  

_cdnaln_protLetters = "ARNDCQEGHILKMFPSTWYV";

_cdnaln_scoreMatrix = pSM2cSM(_cdnaln_protScoreMatrix, _cdnaln_protLetters);

_cdnaln_alnopts = {};
_cdnaln_alnopts ["SEQ_ALIGN_SCORE_MATRIX"] = _cdnaln_scoreMatrix;
_cdnaln_alnopts ["SEQ_ALIGN_GAP_OPEN"] = 40;
_cdnaln_alnopts ["SEQ_ALIGN_AFFINE"] = 1;
_cdnaln_alnopts ["SEQ_ALIGN_GAP_OPEN2"] = 20;
_cdnaln_alnopts ["SEQ_ALIGN_GAP_EXTEND2"] = 1;
_cdnaln_alnopts ["SEQ_ALIGN_GAP_EXTEND"] = 10;
_cdnaln_alnopts ["SEQ_ALIGN_FRAMESHIFT"] = -2*Min(_cdnaln_protScoreMatrix,0);
_cdnaln_alnopts ["SEQ_ALIGN_CODON_ALIGN"] = 1;
_cdnaln_alnopts ["SEQ_ALIGN_NO_TP"] = 1;
_cdnaln_alnopts ["SEQ_ALIGN_CHARACTER_MAP"] = "ACGT";

_cdnaln_partialScoreMatrices = cSM2partialSMs(_cdnaln_scoreMatrix);

_cdnaln_alnopts ["SEQ_ALIGN_PARTIAL_3x1_SCORES"] = _cdnaln_partialScoreMatrices["3x1"];
_cdnaln_alnopts ["SEQ_ALIGN_PARTIAL_3x2_SCORES"] = _cdnaln_partialScoreMatrices["3x2"];

_cdnaln_outstr = "";
_cdnaln_outstr * 256;
_cdnaln_outstr = "[";
// _cdnaln_scores = { 1, Abs( _cdnaln_seqs ) };

for ( _cdnaln_idx = 0; _cdnaln_idx < _cdnaln_numseqs; _cdnaln_idx += 1 )
{
    // align the sequences and get the score
    _cdnaln_inseqs = {{ _cdnaln_refseq, _cdnaln_seqs[ _cdnaln_idx ] }};
    AlignSequences ( _cdnaln_alnseqs, _cdnaln_inseqs, _cdnaln_alnopts );
    _cdnaln_score = ( _cdnaln_alnseqs[0] )[0] / Abs( ( _cdnaln_alnseqs[0] )[1] );
    if ( _cdnaln_dorevcomp ) {
        // if we are going to check the reverse complement:
        // once again align the sequences ( this time with the reverse complement )
        
        _cdnaln_inseqs_rc = {{ _cdnaln_refseq, RevComp( _cdnaln_seqs[ _cdnaln_idx ] ) }};
        AlignSequences( _cdnaln_alnseqs_rc, _cdnaln_inseqs_rc, _cdnaln_alnopts );
        _cdnaln_score_rc = ( _cdnaln_alnseqs_rc[0] )[0] / Abs( ( _cdnaln_alnseqs_rc[0] )[1] );
        
        if ( _cdnaln_score_rc > _cdnaln_score ) {
            // if the reverse complement score is greater than the regular score, use it instead
            _cdnaln_cleanseqs = CleanAlignment( _cdnaln_alnseqs_rc, 0 );
            _cdnaln_score = _cdnaln_score_rc;
        } else {
            // otherwise just the regular score
            _cdnaln_cleanseqs = CleanAlignment(_cdnaln_alnseqs, 0);
        }
    } else {
        // if we are Not checking the reverse complement, just score the result
        _cdnaln_cleanseqs = CleanAlignment(_cdnaln_alnseqs, 0);
    }

    // trim the sequence back to the true length (unpadded)
    // and modify the alignment score to account for this
    _cdnaln_cleanseq = _cdnaln_cleanseqs[ "seq" ];
    _cdnaln_cleanseq = _cdnaln_cleanseq[ 0 ][ ( _cdnaln_truelen - 1 ) ];
    _cdnaln_score = _cdnaln_score * _cdnaln_scoremod;

    if (_cdnaln_idx > 0)
    {
        _cdnaln_outstr * ",";
    }
    _cdnaln_outstr * ( "[\"" + _cdnaln_cleanseq + "\"," + _cdnaln_score + "," + _cdnaln_cleanseqs["overlap"] + "," + (_cdnaln_score/_cdnaln_cleanseqs["overlap"]-computeExcpectedPerBaseScore()) + "]" );
}

_cdnaln_outstr * "]";
_cdnaln_outstr * 0;

