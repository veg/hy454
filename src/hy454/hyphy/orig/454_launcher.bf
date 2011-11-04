/*
 * This is the main launcher for all 454 analyses files.
*/

DO_PIPELINE = 1;

ExecuteAFile ("../Shared/hiv_1_ref_sequences.ibf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"GrabBag.bf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DBTools.ibf");

timer = Time(1);
RequireVersion ("0.9920060815");

fprintf ( stdout, "Provide a quality filtered 454 file:\n" );
fscanf ( PROMPT_FOR_FILE, "Raw", dataFileString );
_in_FilePath = LAST_FILE_PATH;


ExecuteAFile ( "../Shared/chooseGeneticCode.def" );
_in_GeneticCodeTable = modelType;

fprintf ( stdout, "Minimum read length:\n" );
fscanf ( stdin, "Number", _in_minReadL );

ExecuteAFile ( "../Shared/alignmentScoreMatrices/matrixlist.ibf" );
ChoiceList ( modelIdx,"Choose a score matrix",1,SKIP_NONE, aaModelNames);
_in_scoreMatrix = "" + (modelList [modelIdx])["File"];

fprintf ( stdout, "Minimum coverage:\n" );
fscanf ( stdin, "Number", _in_def_min_coverage );

fprintf ( stdout, "Sliding window width for estimation of nucleotide diversity:\n" );
fscanf ( stdin, "Number", _in_def_window_span );

fprintf ( stdout, "Sliding window stride:\n" );
fscanf ( stdin, "Number", _in_def_stride );

fprintf ( stdout, "Minimum number of reads to be considered a variant:\n" );
fscanf ( stdin, "Number", _in_def_minCopyCount );

fprintf ( stdout, "Nucleotide diversity threshold for the estimation of dual/multi infection (recommended 0.025-0.05):\n" );
fscanf ( stdin, "Number", _in_def_nuc_diversity_threshold );

fprintf ( stdout, "Identify drug resistant sites and compensatory mutations (0/1: No/Yes). Note: Only for rt, pol and prrt:\n" );
fscanf ( stdin, "Number", _in_dodr );

if ( _in_dodr ) {
    fprintf ( stdout, "Minimum Stanford Score to consider a drug resistant residue:\n" );
    fscanf ( stdin, "Number", _in_def_min_drugscore );

    fprintf ( stdout, "Minimum coverage for analysis of drug resistant sites:\n" );
    fscanf ( stdin, "Number", _in_def_min_mdr_coverage );
}

ChoiceList ( gene_choice,"Choose a reference sequence",0,SKIP_NONE,RefSeqNames);
gene_map = { Rows ( RefSeqNames), 1 };
for ( _k = 0; _k < Columns ( gene_choice ); _k = _k + 1 ) {
    gene_map [ gene_choice [ _k ] ] = 1;
}

baseFilePath = _in_FilePath;
fprintf (stdout, "Phase 0: Filtering reads based on the reference alignment\n\n" );

/*PH0: Alignment filtering */
first = 1;
for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
    if ( gene_map [_ak] ) {
        dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
        fprintf ( stdout, "Phase 0: Filtering ", dagene,"\n" );
        if ( first ) {
            first = 0;
        }
        else {
            prevgene = (RefSeqNames[lastgene][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
            daDataFile = baseFilePath + "_uds." + prevgene + ".remaining.fas";
            fscanf ( daDataFile, "Raw", dataFileString );
        }
        resultDB = baseFilePath + "_uds." + dagene + ".cache";
        scoreFileString = "../Shared/alignmentScoreMatrices" + DIRECTORY_SEPARATOR + _in_scoreMatrix;

        _options = {};
        ExecuteCommands ( "_options[\"00\"] = \"" + _in_GeneticCodeTable + "\";" );
        ExecuteCommands ( "_options[\"01\"] = \"" + resultDB + "\";" );
        ExecuteCommands ( "_options[\"02\"] = \"" + _ak + "\";" );
        ExecuteCommands ( "_options[\"03\"] = \"" + scoreFileString + "\";" );
        ExecuteCommands ( "_options[\"04\"] = \"" + _in_minReadL + "\";" );

        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        ExecuteAFile ( "454.bf", _options );
        GLOBAL_FPRINTF_REDIRECT = "";
        lastgene = _ak;
    }
}


/*modify gene_map to only those for which reads were found */
for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
    if ( gene_map [_ak] ) {
        dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
        resultDB = baseFilePath + "_uds." + dagene + ".cache";
        DoSQL ( SQL_OPEN, resultDB, DBID );
        _dbRecordCounter = 0;
        DoSQL ( DBID, "SELECT * FROM AA_ALIGNMENT WHERE COVERAGE>0", "return _CountMatchingRecords(0)");
        totalCoverage = 0 + _dbRecordCounter;
        DoSQL ( SQL_CLOSE, "", DBID );
        if ( totalCoverage == 0 ) {
            gene_map [_ak] = 0;
        }
    }
}


/*PH1: Summary statistics of length, depth, and frequency of majority and minority variants*/

fprintf (stdout, "Phase 1: Estimating summary statistics\n" );
for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
    if ( gene_map [ _ak ] ) {
        dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
        resultDB = baseFilePath + "_uds." + dagene + ".cache";
        csvpath = baseFilePath + "_uds." + dagene + ".posreport.csv";
        _options = {};
        ExecuteCommands ( "_options[\"00\"] = \"" + resultDB + "\";" );
        ExecuteCommands ( "_options[\"01\"] = \"" + _in_def_min_coverage + "\";" );
        ExecuteCommands ( "_options[\"02\"] =\"" + csvpath + "\";" );
        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        ExecuteAFile ( "454_reporter.bf", _options );
        GLOBAL_FPRINTF_REDIRECT = "";
    }
}


/*PH2: Sliding window analysis of nucleotide diversity, tree drawing etc */

fprintf (stdout, "Phase 2: Estimating nucleotide diversity in sliding windows\n" );
for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
    if ( gene_map [_ak] ) {
        dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
        resultDB = baseFilePath + "_uds." + dagene + ".cache";
        _options = {};
        ExecuteCommands ( "_options[\"00\"] = \"" + resultDB + "\";" );
        ExecuteCommands ( "_options[\"01\"] = \"" + _in_def_window_span + "\";" );
        ExecuteCommands ( "_options[\"02\"] = \"" + _in_def_stride + "\";" );
        ExecuteCommands ( "_options[\"03\"] = \"" + _in_def_minCopyCount + "\";" );
        ExecuteCommands ( "_options[\"04\"] = \"" + _in_def_min_coverage + "\";" );
        ExecuteCommands ( "_options[\"05\"] = \"" + _in_def_nuc_diversity_threshold + "\";" );
        _options["06"] = "100";
        _options["07"] = "No";
        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        ExecuteAFile ( "454_sliding_window.wbf", _options );
        GLOBAL_FPRINTF_REDIRECT = "";
    }
}


/*PH3: Estimating the number of mutation rate classes across sites */


fprintf (stdout, "Phase 3: Estimating mutation rate classes\n" );
for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
    if ( gene_map [ _ak ] ) {
        dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
        resultDB = baseFilePath + "_uds." + dagene + ".cache";
        _options = {};
        ExecuteCommands ( "_options[\"00\"] = \"" + resultDB + "\";" );
        ExecuteCommands ( "_options[\"01\"] = \"" + _in_def_min_coverage + "\";" );
        ExecuteCommands ( "_options[\"02\"] = \"" + _ak + "\";" );
        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        ExecuteAFile ( "454_variants.bf", _options );
        GLOBAL_FPRINTF_REDIRECT = "";

        _options = {};
        ExecuteCommands ( "_options[\"00\"] = \"" + resultDB + "\";" );
        ExecuteCommands ( "_options[\"01\"] = \"0\";" );
        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        ExecuteAFile ( "454_rateClass_NEB.bf", _options );
        GLOBAL_FPRINTF_REDIRECT = "";
    }
}


/*PH4: Sitewise diversifying/purifying selection analysis */

fprintf (stdout, "Phase 4: Estimating diversifying/purifying selection at sites\n" );
for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
    if ( gene_map [ _ak ] ) {
        dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
        resultDB = baseFilePath + "_uds." + dagene + ".cache";
        _options = {};
        ExecuteCommands ( "_options[\"00\"] = \"" + resultDB + "\";" );
        ExecuteCommands ( "_options[\"01\"] = \"" + _ak + "\";" );
        ExecuteCommands ( "_options[\"02\"] = \"" + _in_def_min_coverage + "\";" );
        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        ExecuteAFile ( "454_FEL.bf", _options );
        GLOBAL_FPRINTF_REDIRECT = "";
    }
}



/*PH5: Drug resistant mutation screening if rt, pr or pol */
if ( _in_dodr ) {
    fprintf (stdout, "Phase 5: Identifying drug resistant variants\n" );
    for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
        if ( gene_map [ _ak ] && ( _ak == 6 || _ak == 7 || _ak > 10 ) ) { /* ie: pr, rt, prrt, prrt, pol : see ../Shared/hiv_1_ref_sequences.ibf */
            dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
            resultDB = baseFilePath + "_uds." + dagene + ".cache";
            _options = {};
            ExecuteCommands ( "_options[\"00\"] = \"" + resultDB + "\";" );
            ExecuteCommands ( "_options[\"01\"] = \"" + _in_def_min_drugscore + "\";" );
            ExecuteCommands ( "_options[\"02\"] = \"" + _in_def_min_mdr_coverage + "\";" );
            ExecuteCommands ( "_options[\"03\"] = \"" + _ak + "\";" );
            GLOBAL_FPRINTF_REDIRECT = "/dev/null";
            ExecuteAFile ("454_MDR_variants.bf", _options );
            GLOBAL_FPRINTF_REDIRECT = "";

            _options = {};
            ExecuteCommands ( "_options[\"00\"] = \"" + resultDB + "\";" );
            ExecuteCommands ( "_options[\"01\"] = \"1\";" );
            GLOBAL_FPRINTF_REDIRECT = "/dev/null";
            ExecuteAFile ( "454_rateClass_NEB.bf", _options );
            GLOBAL_FPRINTF_REDIRECT = "";


        }
    }
}

/*PH6: Compensatory mutation analysis for (N)NRTI's */
if ( _in_dodr ) {
    fprintf (stdout, "Phase 6: Identifying drug accessory mutations\n" );
    for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
        if ( gene_map [ _ak ] && ( _ak == 7 || _ak > 10 ) ) { /*only rt, prrt, prrt, pol: see ../Shared/hiv_1_ref_sequences.ibf */
            dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
            resultDB = baseFilePath + "_uds." + dagene + ".cache";
            _options = {};
            ExecuteCommands ( "_options[\"00\"] = \"" + _in_GeneticCodeTable + "\";" );
            ExecuteCommands ( "_options[\"01\"] = \"" + _in_scoreMatrix + "\";" );
            ExecuteCommands ( "_options[\"02\"] = \"" + _ak + "\";" );
            ExecuteCommands ( "_options[\"03\"] = \"" + resultDB + "\";" );
            GLOBAL_FPRINTF_REDIRECT = "/dev/null";
            ExecuteAFile ( "454_compensatoryMutations.bf", _options );
            GLOBAL_FPRINTF_REDIRECT = "";
        }
    }
}

/*create a legend table*/
for ( _ak = 0; _ak < Rows ( gene_map ); _ak = _ak + 1 ) {
    if ( gene_map [ _ak ] ) {
        dagene = (RefSeqNames[_ak][0]^{{"HXB2_"}{""}})^{{"NL4_3"}{""}};
        resultDB = baseFilePath + "_uds." + dagene + ".cache";
        _options = {};
        ExecuteCommands ( "_options[\"00\"] = \"" + _in_GeneticCodeTable + "\";" );
        ExecuteCommands ( "_options[\"03\"] = \"" + resultDB + "\";" );
        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        ExecuteAFile ( "454_annotate.bf", _options );
        GLOBAL_FPRINTF_REDIRECT = "";
    }
}


