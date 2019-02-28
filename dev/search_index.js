var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Pseudoseq-1",
    "page": "Home",
    "title": "Pseudoseq",
    "category": "section",
    "text": "(Image: Latest Release) (Image: MIT license) (Image: Stable) (Image: Latest) (Image: Pkg Status)"
},

{
    "location": "#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "The Pseudoseq package allows you to simulate DNA sequencing experiments. DNA sequencing experiments are modelled conceptually as a sampling process."
},

{
    "location": "#Install-1",
    "page": "Home",
    "title": "Install",
    "category": "section",
    "text": "You can install Pseudoseq from the julia REPL:using Pkg\nadd(\"https://github.com/bioinfologics/Pseudoseq.jl.git\")"
},

{
    "location": "#Testing-1",
    "page": "Home",
    "title": "Testing",
    "category": "section",
    "text": "Pseudoseq is tested against Julia 1.X on Linux, OS X, and Windows."
},

{
    "location": "man/introduction/#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "man/introduction/#Introduction-to-Pseudoseq-1",
    "page": "Introduction",
    "title": "Introduction to Pseudoseq",
    "category": "section",
    "text": "The Pseudoseq package allows you to simulate DNA sequencing experiments. DNA sequencing experiments are modelled conceptually as a sampling process."
},

{
    "location": "man/introduction/#Basic-workflow-1",
    "page": "Introduction",
    "title": "Basic workflow",
    "category": "section",
    "text": "First, you create a universe of DNA molecules, that is a set of DNA molecules that represents all the DNA molecule fragments present in a sample of extracted DNA.You then create a subsample of DNA fragments from this starting universe, that has the characteristics you desire, by using a set of transformations: You can tag, fragment, and subsample the DNA fragments in a universe, creating a new universe. Preparing the DNA sample for sequencing then can be considered a set of transformations of a starting universe.You then generate a set of reads from the sample, using a paired-end sequencer, a long read sequencer, or a linked read sequencer. You can create a set of reads from single ends of the fragments, paired-end reads, and paired-end reads that are tagged based on their linkage via some long progenitor DNA fragment. After you\'ve created a set of reads, you can mark which bases in your set of reads will be errors.Then finally, you can use your set of reads to generate either an interleaved FASTQ file, or two FASTQ files (one for R1 reads, and one for R2 reads).This manual will present several walkthroughs, showing how to emulate various sequencing experiments using different technologies."
},

{
    "location": "man/pe_walkthrough/#",
    "page": "Paired end reads",
    "title": "Paired end reads",
    "category": "page",
    "text": ""
},

{
    "location": "man/pe_walkthrough/#Walkthrough:-paired-end-reads-1",
    "page": "Paired end reads",
    "title": "Walkthrough: paired end reads",
    "category": "section",
    "text": ""
},

{
    "location": "man/pe_walkthrough/#Creating-the-universe-1",
    "page": "Paired end reads",
    "title": "Creating the universe",
    "category": "section",
    "text": "All DNA sequencing experiments begin with a sample of tissue or cells. Some hair, some blood and so forth.Such a sample undergoes a DNA extraction preparation in the laboratory, after which the genome exists as number of DNA molecules, suspended in a solution.Many cells are typically used as raw input material for DNA extraction, and so the extracted DNA material contains a great many copies of the genome.In Pseudoseq, we call this extracted genetic material (or rather the abstraction of it) the universe. So named, as it is the totality of all copies of the genome that exist in your simulation, from which all subsequent library prep will be done, and from which all reads will be sequenced.Pseudoseq abstracts DNA sequencing experiments as sampling processes, because this is what they are from a statistical point of view: Just as a quadrat placed at random on the forest floor provides a small sample of it\'s species composition, so it is that a sequencing read provides a small sample of the composition of motifs present in a genome. This sample is the universe from which all samples will be drawn.Starting with a FASTA formatted file containing the genome that you want to simulate sequencing for, we create such a universe with the makeuniverse function as follows:universe = makeuniverse(\"mygenome.fasta\", 5000)Where the second argument is the number of copies of the genome you want to exist in your universe. The example above would create a universe of 5000 copies of the genome in \"mygenome.fasta\".tip: Tip\nYou can think of your universe at this stage as containing many copies of the same genome sequence:   1. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   2. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   3. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n                                ...\n4999. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n5000. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC"
},

{
    "location": "man/pe_walkthrough/#Fragmenting-the-molecules-1",
    "page": "Paired end reads",
    "title": "Fragmenting the molecules",
    "category": "section",
    "text": "So we have created a perfect universe of 5000 genome molecules.If DNA sequencing machines could start at one end of a molecule and sequence all the way to the end with reasonable accuracy and throughput, all we would have to do is select a full DNA molecule from our perfect universe, and run it through a sequencing process to get our genome read. Actually, we\'d run several of the molecules from our perfect universe through the sequencer, and use those  reads to form a consensus that eliminates any mistakes a sequencer machine might make - no technology is completely error proof after all and all detection equipment has a certain rate of error.Sadly, we don\'t have such an amazing sequencer technology that allows the reading of a chromosome in its entirety. Instead, shorter reads are taken of shorter fragments of DNA (even long read technology uses fragments and produces reads, much shorter than the full size of a chromosome molecule).A common step in any DNA sequencing experiment, therefore, is to fragment or shear the DNA molecules that are present in the sample.Therefore, a common step in any sequencing experiment simulation constructed with Pseudoseq is to fragment the DNA molecules in a universe into smaller ones.This is achieved with the fragment method.fragment is one of the universe transformations briefly mentioned in the introduction. Other transformations include subsample and tag.Let\'s use it to cut up the molecules in the universe:cut_universe = fragment(universe, 700)The second value of 700 provided above is the desired expected length of the fragments. The fragment method uses this value, and the length of each molecule it is about to cut, to decide how many breakpoints to scatter across said fragment.tip: Tip\nYou can visualise this process like so:From:   1. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   2. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   3. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n                                ...\n4999. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n5000. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGCTo:   1. CGGACTT GAATAGC CCAAA GGTTTCGACACGA TCACGAC ACATAAAT TGGCGGAC TTGAATAGC\n   2. CGGA CTTGAAT AGCCCAAAG GTTTCGAC ACGATCACGACACAT AAATTGGCGGA CTTGA ATAGC\n   3. CGGACTTGA ATAGCC CAAAGGT TTCGACACGAT CACGACACA TAAATT GGCGGACTT GAATAGC\n                                ...\n4999. CGGAC TTGAATA GCCCAAAGGTTT CGACACGA TCACGACACAT AAATTG GCGGACTTG AATAGC\n5000. CGGACTTGAA TAGCCCA AAGGTTTCGA CACGATCAC GACACA TAAATTGGCGG ACTTGAAT AGC"
},

{
    "location": "man/pe_walkthrough/#Subsampling-the-molecules-1",
    "page": "Paired end reads",
    "title": "Subsampling the molecules",
    "category": "section",
    "text": "DNA sequencing experiments are sampling processes, not every one of the  DNA molecules in a sample will be read by the sequencer. So Pseudoseq provides a function to subsample a Universe of DNA fragments.This is done with the sample function:subsampled_universe = subsample(cut_universe, 463967)Where the second numeric argument (463967) is the number of molecules you want to sample.Each fragment in the universe is as likely to be sampled as any other. Fragments are sampled without replacement."
},

{
    "location": "man/pe_walkthrough/#Determining-the-number-of-fragments-to-sample-1",
    "page": "Paired end reads",
    "title": "Determining the number of fragments to sample",
    "category": "section",
    "text": "But how did I get to the number 463967?Prior to sampling, there is a universe of 33139992 DNA fragments, of on average 700bp in length. There is variation - some fragments are bigger, some are smaller.We can estimate how often a base position in the genome is represented in this universe of ~700bp fragments with a simple formula:C = fracLNGWhere C is the expected coverage, L is the average length of the fragments, and N is the number of fragments, and G is the number of bases in the genome.Given that L = 700, N = 33139992, and G = 4639675:C = frac700 times 331399924639675 = 5000Each base in the genome is expected to be represented 5000 times (give or take) in the universe. This makes sense, because if you recall how we created our universe.....universe = makeuniverse(\"mygenome.fasta\", 5000).....we started with 5000 copies of the genome!Ok, so say I wanted to subsample the universe such that each base position in the genome is represented ~50 times. That is, to achieve 50x coverage, you can determine the number of fragments to subsample from the universe, by reversing the formula:N = fracCGLFilling in the relevant values:N = frac50 times 4639675700 = 331405So if 331405 700bp fragments were sampled from the universe of 33139992 700bp fragments, you can expect each base in the genome to be represented approximately 50 times.But wait! It\'s not quite a simple as that! The answer of 331405 would be correct if we were to just pull those fragments out of the universe, and then sequence a complete read of each 700bp fragment in the subsample in its entirety.However, that\'s not what happens in short read sequencing. Instead, some sub-stretch of the fragment is read by the sequencer, and in the case of paired end sequencing, two sub-stretches of each fragment are read by the sequencer: one from each end of the fragment:Read 1\n------->\nTGAATAGCCCAAAGGTTTCGACA\n|||||||||||||||||||||||\nACTTATCGGGTTTCCAAAGCTGT\n              <--------\n                 Read 2Some of the fragment in the middle does not get read by the sequencer.So how did I get the answer 463967 above?First, I know this simulation is of paired end sequencing, and that the length of each read will be 250bp when we get to that stage of the process. We\'re going to sample fragments of 700bp from the universe, and read 250bp in from either end.I can estimate how many 700bp fragments I need to sample to achieve 50x coverage with a read length of 250bp, by using the formula with a revised L, setting it to the read length 250bp, instead of the fragment length of 700bp:N = frac50 times 4639675250This gives me an N of 927935, which I must divide by 2 (I\'m going to be reading two ends of each 700bp DNA fragment I sample), which gives me 463967.tip: Tip\nPseudoseq provides a helper function that assists in this type of calculation:genome_size = 4639675\nexpected_coverage = 50\nread_length = 250\n\nN = needed_sample_size(expected_coverage, genome_size, readlength)\n\n# Divide by 2 as we\'re doing paired end sequencing i.e 2 reads from one fragment.\ndiv(N, 2)So after all that, we have a subsample of a universe of 700bp fragments. The next step is to produce a set of reads from our universe sample."
},

{
    "location": "man/pe_walkthrough/#Creating-the-reads-1",
    "page": "Paired end reads",
    "title": "Creating the reads",
    "category": "section",
    "text": "We now want to create a set of paired-end reads. We want our reads to be 250bp in length.We use the make_reads function to achieve this:pe_reads = make_reads(PairedEnd, subsampled_universe, 250)The first argument to the function is the read type. Possible read types currently provided are PairedEnd, SingleEnd, and TaggedPairs.Now we have some reads, we should mark positions in the reads that are destined to be errors in the output FASTQ.We do this using the mark_errors function. The function requires an error rate.pe_w_errs = mark_errors(pe_reads, 0.001)"
},

{
    "location": "man/pe_walkthrough/#Generate-FASTQ-files-1",
    "page": "Paired end reads",
    "title": "Generate FASTQ files",
    "category": "section",
    "text": "Now we have some paired end reads and have marked some positions as errors, we can generate FASTQ files. This is done with the generate function.generate(\"myreads.fastq\", pe_w_errs).... aaaand thats it! That\'s how to simulate paired end sequencing and generate read files using Pseudoseq."
},

{
    "location": "man/se_walkthrough/#",
    "page": "Long single end reads",
    "title": "Long single end reads",
    "category": "page",
    "text": ""
},

{
    "location": "man/se_walkthrough/#Walkthrough:-Long,-single-end-reads-1",
    "page": "Long single end reads",
    "title": "Walkthrough: Long, single end reads",
    "category": "section",
    "text": "Let\'s see how you might simulate something like an Oxford Nanopore sequencing experiment."
},

{
    "location": "man/se_walkthrough/#Create-a-universe-1",
    "page": "Long single end reads",
    "title": "Create a universe",
    "category": "section",
    "text": "Let\'s start with a univers just as we did for the paired-end reads walkthrough:universe = makeuniverse(\"mygenome.fasta\", 5000)"
},

{
    "location": "man/se_walkthrough/#Fragmenting-the-molecules-1",
    "page": "Long single end reads",
    "title": "Fragmenting the molecules",
    "category": "section",
    "text": "Now we fragment the molecules in our universe, but this time we will use a much larger fragment size, to simulate higher molecular weight DNA that long read technologies typically require:cut_universe = fragment(universe, 40000)Now we have a universe of large cut fragments. We can do a subsample of the universe. We will estimate the number of fragments to sample using the average length of the molecules (40000). We will use the needed_sample_size function to help us out, as we did in the paired end walkthrough.genome_size = 4639675\nexpected_coverage = 30\nreadlength = 40000\n\nN = needed_sample_size(expected_coverage, genome_size, readlength)\n\nuni_sample = subsample(cut_universe, N)"
},

{
    "location": "man/se_walkthrough/#Creating-the-reads-1",
    "page": "Long single end reads",
    "title": "Creating the reads",
    "category": "section",
    "text": "Ok now we\'ve come to creating the reads, we do this with the make_reads function. By only specifying the read type, and the sample of our universe, and omitting the read length we want. Make reads will generate reads from the entire length of each molecule in our sample. We elected to do this to emulate what nanopore seuqncing is supposed to do: It takes an entire DNA fragment, feeds it through an electrically charged pore, producing a read for the entire fragment. se_reads = make_reads(SingleEnd, uni_sample)Now we have some reads, we should mark positions in the reads that are destined to be errors in the output FASTQ.We do this using the mark_errors function. The function requires an error rate.se_w_errs = mark_errors(se_reads, 0.1)We\'ve marked many more positions as errors than we did for the paired end tutorial, this is because current long read sequencing technology has a much higher error rate and short read sequencing technologies, in general.Ok, now we have marked the errors, we can generate a FASTQ file with these reads:generate(\"longreads.fastq\", se_w_errs)And thats all there is to it!"
},

{
    "location": "man/tg_walkthrough/#",
    "page": "Tagged paired end reads",
    "title": "Tagged paired end reads",
    "category": "page",
    "text": ""
},

{
    "location": "man/tg_walkthrough/#Walkthrough:-tagged-paired-end-reads-1",
    "page": "Tagged paired end reads",
    "title": "Walkthrough: tagged paired end reads",
    "category": "section",
    "text": "Let\'s see how we might simulate something like an 10x sequencing experiment.This experiment is simulated much in the same way as demonstrated in the paired-end sequencing walkthrough, except that there is an extra fragmentation step and a tagging step. Let\'s begin!"
},

{
    "location": "man/tg_walkthrough/#Creating-the-universe-1",
    "page": "Tagged paired end reads",
    "title": "Creating the universe",
    "category": "section",
    "text": "Let\'s create a universe like the one in the paired-end sequencing walkthrough.universe = makeuniverse(\"mygenome.fasta\", 5000)"
},

{
    "location": "man/tg_walkthrough/#Fragmenting-the-molecules-1",
    "page": "Tagged paired end reads",
    "title": "Fragmenting the molecules",
    "category": "section",
    "text": "Ok, this time we are going to fragment the molecules, but first at a larger size:big_cut_universe = fragment(universe, 40000)Now we have a universe of fragments approximately 40,000bp in size."
},

{
    "location": "man/tg_walkthrough/#Tag-the-large-fragments-1",
    "page": "Tagged paired end reads",
    "title": "Tag the large fragments",
    "category": "section",
    "text": "Ok, now we will tag these large fragments randomly, using a pool of N distinct tags. Once you tag a fragment in a universe, any other fragments that are derived from that tagged fragment will inherit the same tag. For example, if a DNA fragment in a universe is tagged, and then it is subsequently fragmented during a fragment transform, then all the smaller fragments derived from that long fragment will inherit that long fragments tag. This allows shorter fragments to possess longer range information in the form of these tags, and this is the basis of 10x and similar technologies.tagged_universe = tag(big_cut_universe, 1000000)Here I\'m going to use a pool of 1,000,000 distinct tags. Which fragment gets a certain tag is random. The size of the tag pool, and the number of fragments in your universe will determine how likely it is that any two fragments get the same tag."
},

{
    "location": "man/tg_walkthrough/#Fragmenting-the-molecules,-again-1",
    "page": "Tagged paired end reads",
    "title": "Fragmenting the molecules, again",
    "category": "section",
    "text": "Ok now this time we will fragment the large tagged molecules in our universe down to an average expected fragment size of 700bp, just as we did when we were simulating paired-end sequencing.cut_universe = fragment(tagged_universe, 700)"
},

{
    "location": "man/tg_walkthrough/#Subsampling-the-molecules-1",
    "page": "Tagged paired end reads",
    "title": "Subsampling the molecules",
    "category": "section",
    "text": "Ok, now let\'s subsample the universe just as we did for the paired end sequencing simulation, since the genome size, and read length is the same as in that simulation. We already know to sample 463967 fragments in order to achieve approximately 50x coverage:subsampled_universe = subsample(cut_universe, 463967)"
},

{
    "location": "man/tg_walkthrough/#Creating-the-reads-1",
    "page": "Tagged paired end reads",
    "title": "Creating the reads",
    "category": "section",
    "text": "We now want to create a set of paired-end reads. We want our reads to be 250bp in length.We use the make_reads function to achieve this:tagged_reads = make_reads(TaggedPairs, subsampled_universe, 250)Now we have some reads, we should mark positions in the reads that are destined to be errors in the output FASTQ.We do this using the mark_errors function. The function requires an error rate.tagged_w_errs = mark_errors(tagged_reads, 0.001)Currently, every position in every read is equally likely to be marked as an error."
},

{
    "location": "man/tg_walkthrough/#Generate-FASTQ-files-1",
    "page": "Tagged paired end reads",
    "title": "Generate FASTQ files",
    "category": "section",
    "text": "Now we have some paired end reads and have marked some positions as errors, we can generate FASTQ files. This is done with the generate function.generate(\"tagged_reads.fastq\", tagged_w_errs)"
},

{
    "location": "api/universe/#",
    "page": "Universe",
    "title": "Universe",
    "category": "page",
    "text": ""
},

{
    "location": "api/universe/#API:-Universe-1",
    "page": "Universe",
    "title": "API: Universe",
    "category": "section",
    "text": ""
},

{
    "location": "api/universe/#Pseudoseq.makeuniverse",
    "page": "Universe",
    "title": "Pseudoseq.makeuniverse",
    "category": "function",
    "text": "makeuniverse(gen::Vector{BioSequence{DNAAlphabet{2}}}, ng::Int = 1, iscircular::Bool = false)\n\nCreate a universe of ng copies of a genome defined by the gen vector of sequences.\n\nnote: Note\nThe argument iscircular is currently not used.\n\n\n\n\n\nmakeuniverse(rdr::FASTA.Reader, ng::Int = 1, iscircular::Bool = false)\n\nCreate a universe of ng copies of the genome read in from the FASTA.Reader.\n\nnote: Note\nThe argument iscircular is currently not used.\n\n\n\n\n\nmakeuniverse(file::String, ng::Int, iscircular::Bool = false)\n\nCreate a universe of ng copies of the genome in the fasta formatted file.\n\nnote: Note\nThe argument iscircular is currently not used.\n\n\n\n\n\n"
},

{
    "location": "api/universe/#Pseudoseq.fragment-Tuple{Pseudoseq.Universe,Int64}",
    "page": "Universe",
    "title": "Pseudoseq.fragment",
    "category": "method",
    "text": "fragment(u::Universe, meansize::Int)\n\nCreate a new universe by breaking up the DNA fragments in an input universe.\n\nThis method breaks up a DNA molecule in a universe u, such that the average length of the fragments is approximately meansize.\n\nIt fragments a molecule by scattering an appropriate number of breakpoints across the molecule, before cutting the molecule at those breakpoints.\n\nnote: Note\nBreakpoints are scattered entirely at random across a molecule. No two or more breakpoints can fall in exactly the same place, as those positions are sampled without replacement.\n\nnote: Note\nThe appropriate number of breakpoints to scatter across a molecule is calculated as:fracLS - 1Where L is the length of the molecule being fragmented, and S is the desired expected fragment size. This calculation assumes breakpoints fall randomly across the molecule (see above note).\n\nnote: Note\nIf a DNA molecule being fragmented is smaller than the desired meansize, then it will not be broken, it will simply be included in the new universe.\n\n\n\n\n\n"
},

{
    "location": "api/universe/#Pseudoseq.subsample-Tuple{Pseudoseq.Universe,Int64}",
    "page": "Universe",
    "title": "Pseudoseq.subsample",
    "category": "method",
    "text": "subsample(u::Universe, n::Int)\n\nCreate a new universe by sampling an input universe.\n\nnote: Note\nDNA molecules in the input universe u are selected according to the uniform distribution; no one molecule is more or less likely to be selected than another.\n\nnote: Note\nSampling is done without replacement, so it is impossible for the new universe that is created to recieve one molecule the input universe twice.\n\n\n\n\n\n"
},

{
    "location": "api/universe/#Exported-1",
    "page": "Universe",
    "title": "Exported",
    "category": "section",
    "text": "makeuniverse\nfragment(u::Pseudoseq.Universe, meansize::Int)\nsubsample(u::Pseudoseq.Universe, n::Int)"
},

{
    "location": "api/reads/#",
    "page": "Reads",
    "title": "Reads",
    "category": "page",
    "text": ""
},

{
    "location": "api/reads/#API:-Reads-1",
    "page": "Reads",
    "title": "API: Reads",
    "category": "section",
    "text": ""
},

{
    "location": "api/reads/#Pseudoseq.make_reads",
    "page": "Reads",
    "title": "Pseudoseq.make_reads",
    "category": "function",
    "text": "make_reads(::Type{PairedEnd}, u::Universe, flen::Int, rlen::Int = flen)\n\nCreate a set of paired-end reads from a universe of DNA molecules u.\n\nflen sets the length of forward read, and rlen sets the length of the reverse read. If you only provide flen, then the function sets rlen = flen.\n\nnote: Note\nIf a molecule in the universe is not long enough to create a forward and/or reverse read, then that molecule will simply be skipped. \n\n\n\n\n\nmake_reads(::Type{SingleEnd}, u::Universe, len::Int)\n\nCreate a set of single-end reads from a universe of DNA molecules u.\n\nlen sets the length of the reads.\n\nThe end (strand) from which the reading begins for each DNA molecule in the universe is determined at random for each molecule, with 50:50 probability.\n\nIf you don\'t provide a value for len, then the function will read each DNA molecule in it\'s entirety.\n\nnote: Note\nIf a molecule in the universe is not long enough to create a forward and/or reverse read, then that molecule will simply be skipped. \n\n\n\n\n\nmake_reads(::Type{TaggedPairs}, u::Universe, flen::Int, rlen::Int = flen)\n\nCreate a set of tagged paired-end reads from a universe of DNA molecules u.\n\nflen sets the length of forward read, and rlen sets the length of the reverse read. If you only provide flen, then the function sets rlen = flen.\n\nnote: Note\nIf a molecule in the universe is not long enough to create a forward and/or reverse read, then that molecule will simply be skipped. \n\n\n\n\n\n"
},

{
    "location": "api/reads/#Pseudoseq.mark_errors",
    "page": "Reads",
    "title": "Pseudoseq.mark_errors",
    "category": "function",
    "text": "mark_errors(reads::Reads, rate::Float64)\n\nCreate a new set of reads, with errors, from an input set of reads.\n\nWhen you first create a set of reads using a make_reads method, all the reads in that set are perfect reads. In real sequencing experiments, sequencers make errors when reading a DNA molecule, and are characterised by an error rate. This function lets you simulate this characteristic of sequencers by marking positions in the reads that are destined to be errors in the output FASTQ.\n\nnote: Note\nCurrently, every position in every read is equally likely to be marked as an error.\n\n\n\n\n\n"
},

{
    "location": "api/reads/#Exported-1",
    "page": "Reads",
    "title": "Exported",
    "category": "section",
    "text": "make_reads\nmark_errors"
},

]}
