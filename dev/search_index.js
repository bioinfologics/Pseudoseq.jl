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
    "text": "Fake genomes, fake sequencing, real insights.(Image: Latest Release) (Image: MIT license) (Image: Stable) (Image: Latest) (Image: Pkg Status) (Image: DOI)"
},

{
    "location": "#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "The Pseudoseq package allows you to build arbitrary genomes, and simulate DNA sequencing experiments. DNA sequencing experiments are modelled conceptually as a sampling process."
},

{
    "location": "#Install-1",
    "page": "Home",
    "title": "Install",
    "category": "section",
    "text": "Pseudoseq is built with BioJulia, and is designed with compatibility with the BioJulia ecosystem of tools in mind. Pseudoseq is made available to install through BioJulia\'s package registry.Julia by default only watches the \"General\" package registry, so before you start, you should add the BioJulia package registry.Start a julia terminal, hit the ] key to enter pkg mode (you should see the prompt change from julia> to pkg>), then enter the following command:registry add https://github.com/BioJulia/BioJuliaRegistry.gitAfter you\'ve added the registry, you can install Pseudoseq from the julia REPL. Press ] to enter pkg mode again, and enter the following:add Pseudoseq"
},

{
    "location": "#Testing-1",
    "page": "Home",
    "title": "Testing",
    "category": "section",
    "text": "Pseudoseq is tested against Julia 1.X on Linux, OS X, and Windows."
},

{
    "location": "man/sequencing/concepts/#",
    "page": "Core concepts & workflow",
    "title": "Core concepts & workflow",
    "category": "page",
    "text": ""
},

{
    "location": "man/sequencing/concepts/#Sequencing:-*Core-concepts-and-basic-workflow*-1",
    "page": "Core concepts & workflow",
    "title": "Sequencing: Core concepts & basic workflow",
    "category": "section",
    "text": "Pseudoseq abstracts DNA sequencing experiments as sampling processes, because this is what they are from a statistical point of view. Just as a quadrat placed at random on the forest floor provides a small sample of it\'s species composition, so it is that a sequencing read provides a small sample of the composition of the motifs present in a genome.This manual includes several examples showing how to emulate various sequencing experiments using different technologies. But the core workflow, and important concepts are outlined below.tip: Tip\nThe user can use Pseudoseq\'s API to script each stage of the flow outlined below themselves, or they can use the sequence function, which is the highest-level user facing function. Every example in this section of the manual, will show you how to use both options to achieve the same goal.Anyway, the core sequencing workflow in Pseudoseq is as follows..."
},

{
    "location": "man/sequencing/concepts/#.-Create-a-pool-of-DNA-molecules-1",
    "page": "Core concepts & workflow",
    "title": "1. Create a pool of DNA molecules",
    "category": "section",
    "text": "In reality, all DNA sequencing experiments begin with a sample of tissue or cells. DNA is extracted from the sample in the laboratory, after the genomes of the cells exist as number of DNA molecules, suspended in a solution.In Pseudoseq, such a collection of DNA molecules is called the molecule pool, it is created with the makepool function.In the beginning of a Psueodseq simulation script you create a pool that is the totality of all copies of the genome that exist in your simulation.tip: Tip\nYou can think of the pool at this stage as containing many copies of the same genome sequence:   1. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   2. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   3. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n                                ...\n4999. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n5000. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC"
},

{
    "location": "man/sequencing/concepts/#.-Process-the-DNA-molecule-pool-1",
    "page": "Core concepts & workflow",
    "title": "2. Process the DNA molecule pool",
    "category": "section",
    "text": "You then subject your starting DNA molecule pool to a series of transformations, until the pool has the properties you want to emulate.Here we will describe the different transformations you can apply to a DNA molecule pool:"
},

{
    "location": "man/sequencing/concepts/#Fragmenting-the-pool-1",
    "page": "Core concepts & workflow",
    "title": "Fragmenting the pool",
    "category": "section",
    "text": "In an ideal world, if DNA sequencing machines could start at one end of a molecule and read the sequence all the way to the end with reasonable accuracy and throughput, you could simulate that by selecting molecules from your pool, and producing a read file. Give or take some inevitable errors (all detection equipment has a rate of error), assembling a genome would be simple.Sadly, we don\'t have such an amazing sequencer technology that allows the reading of a chromosome in its entirety. Instead, shorter reads are taken of short fragments of DNA. Even long read technology uses fragments and produces reads, much shorter than the full size of a chromosome molecule.A common step in any DNA sequencing experiment, therefore, is to fragment or shear the DNA molecules that are present in the pool.This is achieved in Pseudoseq with the fragment function.tip: Tip\nYou can visualise this process like so:From:   1. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   2. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n   3. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n                                ...\n4999. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGC\n5000. CGGACTTGAATAGCCCAAAGGTTTCGACACGATCACGACACATAAATTGGCGGACTTGAATAGCTo:   1. CGGACTT GAATAGC CCAAA GGTTTCGACACGA TCACGAC ACATAAAT TGGCGGAC TTGAATAGC\n   2. CGGA CTTGAAT AGCCCAAAG GTTTCGAC ACGATCACGACACAT AAATTGGCGGA CTTGA ATAGC\n   3. CGGACTTGA ATAGCC CAAAGGT TTCGACACGAT CACGACACA TAAATT GGCGGACTT GAATAGC\n                                ...\n4999. CGGAC TTGAATA GCCCAAAGGTTT CGACACGA TCACGACACAT AAATTG GCGGACTTG AATAGC\n5000. CGGACTTGAA TAGCCCA AAGGTTTCGA CACGATCAC GACACA TAAATTGGCGG ACTTGAAT AGC"
},

{
    "location": "man/sequencing/concepts/#Subsampling-molecules-from-a-pool-1",
    "page": "Core concepts & workflow",
    "title": "Subsampling molecules from a pool",
    "category": "section",
    "text": "DNA sequencing experiments are sampling processes, not every one of the  DNA molecules in a pool will be read by the sequencer.Randomly subsampling (without replacement) DNA molecules of a pool is achieved in Pseudoseq with the subsample function."
},

{
    "location": "man/sequencing/concepts/#Determining-the-number-of-molecules-to-subsample-1",
    "page": "Core concepts & workflow",
    "title": "Determining the number of molecules to subsample",
    "category": "section",
    "text": "You can estimate how often a base position in the genome is represented in a DNA molecule pool using the following formula.C = fracLNGWhere C is the expected coverage, L is the average length of the fragments, and N is the number of fragments, and G is the number of bases in the genome.If you wanted to subsample a pool such that each base position in a genome is represented ~50 times (i.e. to achieve 50x coverage), you can determine the number of fragments to subsample from the universe, by reversing the formula:N = fracCGLFor example, if you wanted to sequence 250bp paired-end reads, at 50x coverage, with a genome size of 4639675bp:frac50 times 4639675250 = 927935Remembering that you get two reads from one DNA molecule with paired-end sequencing, you know to subsample 927935  2 = 463967 DNA molecules from a pool.tip: Tip\nPseudoseq provides a helper function that assists in this type of calculation:genome_size = 4639675\nexpected_coverage = 50\nread_length = 250\n\nN = needed_sample_size(expected_coverage, genome_size, readlength)\n\n# Divide by 2 as we\'re doing paired end sequencing.\ndiv(N, 2)"
},

{
    "location": "man/sequencing/concepts/#Tagging-molecules-in-a-pool-1",
    "page": "Core concepts & workflow",
    "title": "Tagging molecules in a pool",
    "category": "section",
    "text": "Some DNA sequencing technologies work by allowing short reads to contain longer range information by tagging molecules in a pool. The basic idea being that if two short reads come from the same large DNA molecule, they will have the same tag.If a DNA fragment in a universe is tagged, and then it is subsequently fragmented during a fragment transform, then all the smaller fragments derived from that long fragment will inherit that long fragment\'s tag. This allows shorter fragments to possess longer range information in the form of these tags, and this is the basis of 10x sequencing and similar technologies.Pseudoseq lets you attach tags to molecules in a pool using the tag function."
},

{
    "location": "man/sequencing/concepts/#.-Generating-reads-1",
    "page": "Core concepts & workflow",
    "title": "3. Generating reads",
    "category": "section",
    "text": "Next, you generate a set of reads from your transformed pool.  Pseudoseq allows you to create paired-end sequencing reads, single-end sequencing reads, and linked paired-end sequencing reads.You create a set of reads using the make_reads function.All sequencers have a certain error rate, and so after you\'ve created a set of reads, you can use the mark_errors function to randomly mark a set of bases in your reads that are destined to be errors in the output file.Finally, you use your set of reads to generate either an interleaved FASTQ file, or two FASTQ files (one for R1 reads, and one for R2 reads)."
},

{
    "location": "man/sequencing/examples/pe-example/#",
    "page": "Paired end reads",
    "title": "Paired end reads",
    "category": "page",
    "text": "EditURL = \"https://github.com/bioinfologics/Pseudoseq.jl/blob/master/examples/sequencing/pe-example.jl\""
},

{
    "location": "man/sequencing/examples/pe-example/#Example:-paired-end-sequencing-1",
    "page": "Paired end reads",
    "title": "Example: paired-end sequencing",
    "category": "section",
    "text": "This is an example generated from this source file: pe-example.jl You are seeing the online documentation version. The corresponding notebook can be found here: pe-example.ipynb and the script can be found here: pe-example.jlFor the simulation we are going to:Create a pool of 5000 copies of a reference genome.\nFragment the DNA molecules in the pool, to an average length of 700bp.\nSubsample the molecules in the pool to achieve approximatly 50x coverage.\nCreate a set of 250bp paired-end reads.\nApply errors to the paired-end reads at a rate of 0.001 (.1%).\nGenerate an output FASTQ file.using Pseudoseq"
},

{
    "location": "man/sequencing/examples/pe-example/#Using-the-[sequence](@ref)-method-1",
    "page": "Paired end reads",
    "title": "Using the sequence method",
    "category": "section",
    "text": "First, let\'s see how we do this with the sequence method. The first two parameters we give to the function will be the input genome we want to sequence, and the destination FASTQ file for output reads. Here we are setting:The number of genome copies in the molecule pool to 5000.\nThe average fragment size to 700bp.\nThe sampling coverage to 50x.\nThe read length to 250bp.\nThe per base read error rate to 0.001.\nThe fact we want paired-ends of fragments to be read (paired) to true.sequence(\"ecoli-ref.fasta\", \"pe-reads.fastq\"; ng = 5000, flen = 700, cov = 50, paired = true, rdlen = 250, err = 0.001)"
},

{
    "location": "man/sequencing/examples/pe-example/#Using-the-Pseudoseq-API-1",
    "page": "Paired end reads",
    "title": "Using the Pseudoseq API",
    "category": "section",
    "text": "Here\'s how to achieve the same thing, using the Pseudoseq API. It is nessecery to use the API if you want to do something that is not anticipated by the available functionality of the sequence method: the cost of conveinience is fewer options.Starting with a FASTA formatted file containing the genome we want to sequence, we create a pool with 5000 copies of the genome.pool = makepool(\"ecoli-ref.fasta\", 5000)Next we use the fragment function to make a pool of shorter DNA molecules.cutpool = fragment(pool, 700)We need to determine the number of molecules to sample, and subsample the pool:genome_size = 4639675\nexpected_coverage = 50\nread_length = 250\n\nN = needed_sample_size(expected_coverage, genome_size, read_length)\nN = div(N, 2) # Divide by 2 as we\'re doing paired end sequencing.\n\nsampledpool = subsample(cutpool, N)We now want to create a set of paired-end reads. We want our reads to be 250bp in length.pe_reads = make_reads(PairedEnd, sampledpool, 250)Now we have some reads, we should mark positions in the reads that are destined to be errors in the output FASTQ.pe_w_errs = mark_errors(pe_reads, 0.001)Now we have some paired end reads and have marked some positions as errors, we can generate FASTQ files.generate(\"pe-reads.fastq\", pe_w_errs)#-This page was generated using Literate.jl."
},

{
    "location": "man/sequencing/examples/se-example/#",
    "page": "Long single end reads",
    "title": "Long single end reads",
    "category": "page",
    "text": "EditURL = \"https://github.com/bioinfologics/Pseudoseq.jl/blob/master/examples/sequencing/se-example.jl\""
},

{
    "location": "man/sequencing/examples/se-example/#Example:-long,-single-end-reads-1",
    "page": "Long single end reads",
    "title": "Example: long, single end reads",
    "category": "section",
    "text": "This is an example generated from this source file: se-example.jl You are seeing the online documentation version. The corresponding notebook can be found here: se-example.ipynb and the script can be found here: se-example.jlLet\'s see how you might simulate something like an Oxford Nanopore sequencing experiment.For the simulation we are going to:Create a pool of 5000 copies of a reference genome.\nFragment the DNA molecules in the pool, to an average length of 40,000bp.\nSubsample the molecules in the pool to achieve approximatly 30x coverage.\nCreate a set of single-end reads, the enitre length of each molecule.\nApply errors to the reads at a rate of 0.10 (1 error every 10bp).\nGenerate an output FASTQ file.using Pseudoseq"
},

{
    "location": "man/sequencing/examples/se-example/#Using-the-[sequence](@ref)-method-1",
    "page": "Long single end reads",
    "title": "Using the sequence method",
    "category": "section",
    "text": "First, let\'s see how we do this with the sequence method. The first two parameters we give to the function will be the input genome we want to sequence, and the destination FASTQ file for output reads. Here we are setting:The number of genome copies in the molecule pool to 5000.\nThe average fragment size to 40000bp.\nThe sampling coverage to 30x.\nThe read length to nothing, which will make the sequencer read the whole length of any DNA fragment.\nThe per base read error rate to 0.1.\nThe fact we want paired-ends of fragments to be read (paired) to false.sequence(\"ecoli-ref.fasta\", \"longreads.fastq\"; ng = 5000, flen = 40000, cov = 30, rdlen = nothing, err = 0.1, paired = false)"
},

{
    "location": "man/sequencing/examples/se-example/#Using-the-Pseudoseq-API-1",
    "page": "Long single end reads",
    "title": "Using the Pseudoseq API",
    "category": "section",
    "text": "Here\'s how to achieve the same thing, using the Pseudoseq API. It is nessecery to use the API if you want to do something that is not anticipated by the available functionality of the sequence method: the cost of conveinience is fewer options.Let\'s start with a pool of 5000 copies of a genome contained in a FASTA file:pool = makepool(\"ecoli-ref.fasta\", 5000)Cut the pool of DNA into fragments of an average length of 40,000bpcutpool = fragment(pool, 40000)Now we\'ll estimate the number of fragments we need to sample from the pool to achieve 30x coverage.genome_size = 4639675\nexpected_coverage = 30\nreadlength = 40000\n\nN = needed_sample_size(expected_coverage, genome_size, readlength)\n\nsampledpool = subsample(cutpool, N)By using the make_reads function without specifying a read length, the function will generate reads from the entire length of each molecule in the pool. We do this to emulate what Nanopore sequencing is supposed to do: It takes an entire DNA fragment, feeds it through an electrically charged pore, producing a read for the entire fragment.se_reads = make_reads(SingleEnd, sampledpool)Long read sequencer have much higher error rates than short read sequencers so we use a error rate of 0.1.se_w_errs = mark_errors(se_reads, 0.1)Finally produce the ouput FASTQ file.generate(\"longreads.fastq\", se_w_errs)#-This page was generated using Literate.jl."
},

{
    "location": "man/sequencing/examples/tg-example/#",
    "page": "Tagged paired end reads",
    "title": "Tagged paired end reads",
    "category": "page",
    "text": "EditURL = \"https://github.com/bioinfologics/Pseudoseq.jl/blob/master/examples/sequencing/tg-example.jl\""
},

{
    "location": "man/sequencing/examples/tg-example/#Example:-tagged-paired-end-reads-1",
    "page": "Tagged paired end reads",
    "title": "Example: tagged paired-end reads",
    "category": "section",
    "text": "This is an example generated from this source file: tg-example.jl You are seeing the online documentation version. The corresponding notebook can be found here: tg-example.ipynb and the script can be found here: tg-example.jlLet\'s see how we might simulate something like an 10x sequencing experiment.For this simulation script we will:Create a pool of 5000 copies of a reference genome.\nFragment the DNA molecules in the pool, to an average length of 40,000bp.\nTag the long molecules in the pool randomly with a set of 1,000,000 tags.\nFragment the molecules in the pool to an average length of 700bp.\nSubsample the molecules in the pool to achieve approximatly 50x coverage.\nCreate a set of 250bp paired-end reads.\nApply errors to the paired-end reads at a rate of 0.001 (.1%).\nGenerate an output FASTQ file.using Pseudoseq"
},

{
    "location": "man/sequencing/examples/tg-example/#Using-the-[sequence](@ref)-method-1",
    "page": "Tagged paired end reads",
    "title": "Using the sequence method",
    "category": "section",
    "text": "First, let\'s see how we do this with the sequence method. The first two parameters we give to the function will be the input genome we want to sequence, and the destination FASTQ file for output reads. Here we are setting:The number of genome copies in the molecule pool to 5,000.\nThe number of possible tags to one million.\nThe average fragment size, prior to tagging, to 40,000bp.\nThe average fragment size after tagging, to 700bp.\nThe sampling coverage to 50x.\nThe read length to 250bp.\nThe per base read error rate to 0.001.\nThe fact we want paired-ends of fragments to be read (paired) to true, which is the default.sequence(\"ecoli-ref.fasta\", \"tagged_reads.fastq\"; ng = 5000, tusize = 1000000, taggedflen = 40000, flen = 700, cov = 50, rdlen = 250, err = 0.1)"
},

{
    "location": "man/sequencing/examples/tg-example/#Using-the-Pseudoseq-API-1",
    "page": "Tagged paired end reads",
    "title": "Using the Pseudoseq API",
    "category": "section",
    "text": "Here\'s how to achieve the same thing, using the Pseudoseq API. It is nessecery to use the API if you want to do something that is not anticipated by the available functionality of the sequence method: the cost of conveinience is fewer options.Let\'s start with a pool of 5000 copies of a genome contained in a FASTA file:dnapool = makepool(\"ecoli-ref.fasta\", 5000)Now let\'s cut up the molecules to an average length of 40,000bpcutpool = fragment(dnapool, 40000)Ok, now we will tag these large fragments randomly. Once you tag a fragment in a universe, any other fragments that are derived from that tagged fragment will inherit the same tag.taggedpool = tag(cutpool, 1000000)Here I\'m going to use a pool of 1,000,000 distinct tags. Which fragment gets a certain tag is random. The size of the tag pool, and the number of fragments in your universe will determine how likely it is that any two fragments get the same tag. Now we\'ll fragment the pool againtaggedcutpool = fragment(taggedpool, 700)Subsample the pool of tagged molecules.genome_size = 4639675\nexpected_coverage = 50\nread_length = 250\n\nN = needed_sample_size(expected_coverage, genome_size, read_length)\nN = div(N, 2) # Divide by 2 as we\'re doing paired end sequencing.\n\nsampledpool = subsample(taggedcutpool, N)Now let\'s make some 250bp tagged paired reads and generate some erroneous positions.tagged_reads = make_reads(TaggedPairs, sampledpool, 250)\ntagged_w_errs = mark_errors(tagged_reads, 0.001)Output to FASTQ:generate(\"tagged_reads.fastq\", tagged_w_errs)#-This page was generated using Literate.jl."
},

{
    "location": "man/build-a-genome/concepts/#",
    "page": "Core concepts & workflow",
    "title": "Core concepts & workflow",
    "category": "page",
    "text": ""
},

{
    "location": "man/build-a-genome/concepts/#Build-a-Genome:-*Core-concepts-and-basic-workflow*-1",
    "page": "Core concepts & workflow",
    "title": "Build-a-Genome: Core concepts & basic workflow",
    "category": "section",
    "text": "Pseudoseq allows you to plan and build genomes and chromosomes that have a certain set of features and peculiarities. The purpose for doing this is not to recreate biology perfectly. The purpose is to create genomes you understand fully (where the repeated content is, which positions are heterozygous and so on).Using such genomes can help you both understand and develop an intuition of what current genome assembly tools are doing, and also to help design assembly tools, and perhaps even plan sequencing experiments and form hypotheses.This manual includes several examples showing how to plan genomes with certain characteristics. But the core workflow, and important concepts are explained below, in 3 steps."
},

{
    "location": "man/build-a-genome/concepts/#.-Creating-chromosome-blueprints-1",
    "page": "Core concepts & workflow",
    "title": "1. Creating chromosome blueprints",
    "category": "section",
    "text": "Chromosome blueprints are the backbone of simulating genomes with Pseudoseq.Chromosome blueprints determine the nature of one chromosome in a genome.You can think of chromosome blueprints in Pseudoseq as: a collection of operations which, when applied to some seed sequence, result in a set of N homologous sequences.note: Note\nChromosome blueprints are immutable: Adding an operation to a chromosome blueprint creates a new blueprint, which is a copy of the old blueprint with the addition of the new operation.note: Note\nIdeally, for any chromosome blueprint you construct with Pseudoseq, for each operation is possesses, it must be possible (at least in principle) to be able to intuit what the effect of that operation will be on:The Khmer distribution produced by sequencing reads of the fabricated chromosome.\nThe structure of a sequence graph produced by sequencing reads of the fabricated chromosome.Therefore, we have made the design decision that no two operations in a chromosome blueprint may affect the same position(s) of the genome in a conflicting manner. To meet this requirement, certain operations \"consume\" a region of the chromosome planned in a blueprint. If a region is consumed, another operation that would also affect that region cannot be added to the blueprint.Depending on the genome, any given chromosome may be present in multiple copies. Diploids, for example have two copies of every chromosome.The first step in simulating any artificial genome is to create one or more blank chromosome blueprints. The plan_chrom function is used for this.For example, this:using Pseudoseq\nc = plan_chrom(100, 2)will create a blank blueprint for 2 copies of a chromosome of 100bp length. From the output you can see that it prints for you the number of chromosome copies, the length of the chromosomes, and a list of available, unconsumed regions of the chromosome (see note above)."
},

{
    "location": "man/build-a-genome/concepts/#.-Adding-planned-features-to-a-chromosome-blueprint-1",
    "page": "Core concepts & workflow",
    "title": "2. Adding planned features to a chromosome blueprint",
    "category": "section",
    "text": "After creating a fresh chromosome blueprint, no plans (operations) have been added yet.If you were to fabricate this blank blueprint, you would get N identical DNA sequences as output, where N is the number of copies the blueprint was planning.Once you have one or more chromosome blueprints, you can add features to them.This is done with a series of consistently named plan_* functions.note: Note\nRemember; blueprints are immutable, so every time one of these plan_* functions is used to add a feature to a chromosome blueprint, a new chromosome blueprint is created."
},

{
    "location": "man/build-a-genome/concepts/#Repetitions-1",
    "page": "Core concepts & workflow",
    "title": "Repetitions",
    "category": "section",
    "text": "A repetition is a segment of a sequence, that has the exact same motif, as  another segment of the sequence.In Pseudoseq, to plan a repetition, you specify a region the repetition will copy, sometimes called the from region. You also specify a region where the motif in the from region will be replicated, called the to region. So you might find it helpful to imagine a planned repetition as a kind of copy-paste operation that occurs during fabricate.note: Note\nRepetitions consume the to region of the chromosome blueprint to which they are added. Repetitions do not consume the from region, so other operations are free to affect the motif in the from region. Just remember that the repetition will replicate anything in the from region, including other features such as heterozygosity that occur in from.Repetitions are added to a genome blueprint using the plan_repetition function."
},

{
    "location": "man/build-a-genome/concepts/#Heterozygosity-1",
    "page": "Core concepts & workflow",
    "title": "Heterozygosity",
    "category": "section",
    "text": "A heterozygosity describes a base position at which the copies of the chromosome in the blueprint differ.For a blueprint with 2 copies, both copies will differ at a given position.For a triploid, at a heterozygous position, all 3 copies might differ from each other. Alternatively, it is possible that 2 copies are the same, but they differ from a 3rd copy. This applies for blueprints with higher copy numbers too: at a heterozygous position some copies will differ from each other at that position, but some of the copies might be the same.note: Note\nHeterozygosity operations consume the position at which they are defined.Heterozygous positions are planned using the plan_het function.For example, below will make it so as about 50% of the two chromosome copies  are heterozygous:using Pseudoseq\nc = plan_chrom(100, 2)\nchet = plan_het(c, .50, 2)The above use of plan_het allows the function to choose which sites in the blueprint are heterozygous, and how to allocate the bases at the heterozygous sites. We only instruct the function that 50% of the genome should be heterozygous, and that there should be 2 possible bases at each position (the only option really: the blueprint plans for a diploid - 2 chromosome copies).You can also take more fine-grained control:using Pseudoseq, BioSequences\nc = plan_chrom(100, 2)\nchet = plan_het(c, 20, [DNA_T, DNA_A])Here I planned a heterozygous position, at one site in the chromosome (20), and I set the state of the first copy to DNA_T, and the state of the second copy to DNA_A.So as you can see, plan_het is very flexible. Check it\'s API documentation, and the \"Build-A-Yeast\" example walkthrough to see more examples of plan_het use."
},

{
    "location": "man/build-a-genome/concepts/#Utility-functions-1",
    "page": "Core concepts & workflow",
    "title": "Utility functions",
    "category": "section",
    "text": "There is a utility function suggest_regions available to help you plan where to place features.Say you wanted to see where you could place 3 5bp repetitions, you could do the following:r = suggest_regions(chet, 5, 6)So now you have 6 5bp regions, every 2 regions defining a from and a to region for a single repetition. You can provide r as an input to plan_repetition.crep = plan_repetition(chet, r)Another utility function suggest_alleles is available to help you plan nucleotide patterns at heterozygous sites.Say you wanted to plan a pattern in which two of three chromosome copies had the same base, and a third differed.suggest_alleles(3, 2)The above asks for an allele pattern for 3 copies of a chromosome, with 2 possible alleles, which chromosome copy gets which allele is random.You can also take more control and specify which chromosome copy gets which state:suggest_alleles(chet, [1, 1, 2])\n# Or alternatively put:\nsuggest_alleles(chet, [1, 2], [3])In the above, which chromosomes get the same allele is user-specified, but which bases those alleles are, is randomly determined.You can get many suggestions at once:suggest_alleles(5, [1, 1, 2])See the API docs for suggest_alleles for more details on the arguments permitted by the different methods. "
},

{
    "location": "man/build-a-genome/concepts/#.-Fabricate-a-FASTA-from-chromosome-blueprints-1",
    "page": "Core concepts & workflow",
    "title": "3. Fabricate a FASTA from chromosome blueprints",
    "category": "section",
    "text": "Once you have a set of chromosome blueprints with the features planned that you desire, you can fabricate the sequences of these chromosomes.To do that, use the fabricate method.You can simply fabricate the sequences, for use in the interactive session:fabricate(chet)Or have the sequences output to a FASTA formatted file.fabricate(\"mychrom.fasta\", chet)A randomly generated seed sequence is used to start the fabrication process unless you provide one. See the API docs for fabricate for more details.And that\'s all there is to building a chromosome with Pseudoseq. To build a genome with more than one chromosome, simply build a set of chromosome blueprints.Now check out the examples to see how various genomes can be built with Pseudoseq."
},

{
    "location": "man/build-a-genome/examples/build-a-yeast/#",
    "page": "Build-a-Yeast",
    "title": "Build-a-Yeast",
    "category": "page",
    "text": "EditURL = \"https://github.com/bioinfologics/Pseudoseq.jl/blob/master/examples/build-a-genome/build-a-yeast.jl\""
},

{
    "location": "man/build-a-genome/examples/build-a-yeast/#Example:-Build-a-Yeast-1",
    "page": "Build-a-Yeast",
    "title": "Example: Build-a-Yeast",
    "category": "section",
    "text": "This is an example generated from this source file: build-a-yeast.jl You are seeing the online documentation version. The corresponding notebook can be found here: build-a-yeast.ipynb and the script can be found here: build-a-yeast.jlWe are going to use Pseudoseq to create a simple set of fake genomes, based on chromosome 1 of the reference genome for yeast.Load the sequence of chromosome 1 of the yeast reference from FASTA file.\nCreate blank chromosome blueprints for different genome designs.\nEdit the blueprints to achive a: a. Diploid chromosome b. Triploid chromosome c. Tetraploid chromosomeusing BioSequences, Pseudoseq"
},

{
    "location": "man/build-a-genome/examples/build-a-yeast/#Load-the-seed-sequence-1",
    "page": "Build-a-Yeast",
    "title": "Load the seed sequence",
    "category": "section",
    "text": "First load the sequence from FASTA file. This uses tools from BioSequences, which is a dependency of Pseudoseq and so the user should have it available, as julia\'s package manager would have installed any dependencies. The lines below open a fasta file, with a fasta reader, and load a single FASTA record, and get its sequence.refseq = open(FASTA.Reader, \"yeast-chr1.fasta\") do rdr\n    FASTA.sequence(BioSequence{DNAAlphabet{2}}, read(rdr))\nend\n\nreflen = length(refseq)This results in a 230,218nt long sequence called refseq."
},

{
    "location": "man/build-a-genome/examples/build-a-yeast/#Designing-a-diploid-1",
    "page": "Build-a-Yeast",
    "title": "Designing a diploid",
    "category": "section",
    "text": "First lets make a plan for a diploid genome with 1% heterozygosity. Pseudoseq has a lot of flexibilty in how you can define the heterozygosity, but for ease here I\'m going to let it decide where the sites should be, and what the alleles at each site should be, by calling plan_het, providing just the proportion of sites I want to be heterozygous, and the number of alleles I want at each of these heterozygous sites (2).dploid = plan_chrom(reflen, 2)\ndploidhet = plan_het(dploid, .01, 2)Now fabricate the FASTA, with fabricate, providing a filename, and then a series of blueprint => seed pairs. In this case we just give 1 blueprint-seed# pair.\n\nfabricate(\"build-a-yeast-di.fasta\", dploidhet => refseq)Thats it for the diploid genome."
},

{
    "location": "man/build-a-genome/examples/build-a-yeast/#Designing-a-triploid-1",
    "page": "Build-a-Yeast",
    "title": "Designing a triploid",
    "category": "section",
    "text": "Now let\'s make a triploid genome. This should also be simple enough. It\'s much like the process for making the diploid genome, but specifying 3 chromosome copies instead of 2.trploid = plan_chrom(reflen, 3)When I add some heterozygous sites, just like with the diploid case, I\'m going to simply ask for a proportion of heterozygous sites (and let Pseudoseq decide where they should occur. Like last time also I\'m going to make all the sites have 2 alleles, so for every heterozygous site, two of the three chromosome copies willl have the same base, and one of the copies will be different. There are other ways to define some heterozygous sites for example all 3 chromosome copies could have different basestrploidhet = plan_het(trploid, .01, 2)Now fabricate the FASTA, with fabricate, providing a filename, and then a series of blueprint => seed pairs. In this case we just give 1 blueprint-seed# pair.\n\nfabricate(\"build-a-yeast-tri.fasta\", trploidhet => refseq)That\'s it for the triploid genome."
},

{
    "location": "man/build-a-genome/examples/build-a-yeast/#Designing-a-tetraploid-genome-1",
    "page": "Build-a-Yeast",
    "title": "Designing a tetraploid genome",
    "category": "section",
    "text": "Lets make a tetraploid genome, that consists of two diploid genomes, an A genome and a B genome. First I\'m going to create a chromosome blueprint with 4 copies:tetploid = plan_chrom(reflen, 4)We want to plan the heterozygosity in this tetraploid such that the number of heterozygous sites within the A and B genomes is ~1%. But we want ~3% of sites to be heterozygous between the A and B genomes.First we will say copies 1 + 2 in our blueprint are the A genome, and copies 3 + 4 are the B genome. Next consider heterozygous sites with 2 different alleles, A and B, 6 different patterns are possible:Copy 1 Copy 2 Copy 3 Copy 4\nA A B B\nB B A A\nA B A B\nB A B A\nA B B A\nB A A BThe first 2 patterns in the table constitute fixed mutations between genomes A and B, but not polymorphism within genomes A and B. This is because Copy 1 and 3 differ, as do copies 2 and 4, yet copies 1 and 2 are the same, as are copies 3 and 4. To achieve a ~3% divergence between genomes A and B, let\'s plan some divergent mutations by grouping the copy numbers together when planning some heterozygous sites:thet = plan_het(tetploid, .03, [1, 2], [3, 4])By using the groups [1, 2] and [3, 4] we ensure chromosome copy 1 and 2 will have the same allele, as will chromosomes 3 and 4Next we\'ll just add some of the other heterozygosity patterns to achieve the 1% internal heterozygosity desired. We\'ll use a heterozygosity pattern with two states.thet = plan_het(thet, .01, [1, 3], [2, 4])\n\nfabricate(\"build-a-yeast-tet.fasta\", thet => refseq)This page was generated using Literate.jl."
},

{
    "location": "api/sequence/#",
    "page": "sequence",
    "title": "sequence",
    "category": "page",
    "text": ""
},

{
    "location": "api/sequence/#API:-sequence-1",
    "page": "sequence",
    "title": "API: sequence",
    "category": "section",
    "text": ""
},

{
    "location": "api/sequence/#Pseudoseq.sequence",
    "page": "sequence",
    "title": "Pseudoseq.sequence",
    "category": "function",
    "text": "sequence(input, output = nothing; kwargs...)\n\nRun a sequencing experiment from start to finish.\n\nThis method is the highest-level function for running sequencing simulations with Pseudoseq. It runs the Pseudoseq sequencin worksflow described in the manual, tuned by the following keyword parameters:\n\nGlobal parameters:\n\ninput: A filename or FASTA.Reader providing the input genome.\noutput: A filename or FASTA.Writer providing a destination for the output reads.\n\nParameters for workflow step 1 (Create a pool of DNA molecules):\n\nng: Integer; the number of genomes the molecule pool is initialized with.\n\nParameters for step 2 (Processing the DNA molecule pool):\n\ntusize: Integer; the number of possible tags a DNA molecule may be tagged with (default: 0).\ntaggedflen: Integer; the desired average fragment length at which molecules are tagged (default: 0bp). \nflen: Integer; the desired average fragment length at which molecules are sampled (default: 700bp).\ncov: Integer; the desired expected coverage (default: 30x coverage).\n\nParameters for step 3 (Generating reads):\n\npaired: true or false; Whether or not to sequence from both end of each molecule (default: true).\nrdlen: Integer; the desired read length (default: 250bp).\nerr: Float; the desired per-base error rate (default: 0.001).\n\n\n\n\n\n"
},

{
    "location": "api/sequence/#Exported-functions-1",
    "page": "sequence",
    "title": "Exported functions",
    "category": "section",
    "text": "sequence"
},

{
    "location": "api/pool/#",
    "page": "Molecule Pool",
    "title": "Molecule Pool",
    "category": "page",
    "text": ""
},

{
    "location": "api/pool/#API:-MoleculePool-1",
    "page": "Molecule Pool",
    "title": "API: MoleculePool",
    "category": "section",
    "text": ""
},

{
    "location": "api/pool/#Exported-functions-1",
    "page": "Molecule Pool",
    "title": "Exported functions",
    "category": "section",
    "text": ""
},

{
    "location": "api/pool/#Pseudoseq.makepool",
    "page": "Molecule Pool",
    "title": "Pseudoseq.makepool",
    "category": "function",
    "text": "makepool(gen::Vector{BioSequence{DNAAlphabet{2}}}, ng::Int = 1)\n\nCreate a pool of ng copies of a genome defined by the gen vector of sequences.\n\n\n\n\n\nmakepool(rdr::FASTA.Reader, ng::Int = 1)\n\nCreate a pool of ng copies of the genome read in from the FASTA.Reader.\n\n\n\n\n\nmakepool(file::String, ng::Int)\n\nCreate a pool of ng copies of the genome in the fasta formatted file.\n\nnote: Note\nThe argument iscircular is currently not used.\n\n\n\n\n\n"
},

{
    "location": "api/pool/#Making-a-pool-of-molecules-1",
    "page": "Molecule Pool",
    "title": "Making a pool of molecules",
    "category": "section",
    "text": "makepool"
},

{
    "location": "api/pool/#Molecule-pool-transformations-1",
    "page": "Molecule Pool",
    "title": "Molecule pool transformations",
    "category": "section",
    "text": ""
},

{
    "location": "api/pool/#Pseudoseq.fragment-Tuple{Pseudoseq.MoleculePool,Int64}",
    "page": "Molecule Pool",
    "title": "Pseudoseq.fragment",
    "category": "method",
    "text": "fragment(p::MoleculePool, meansize::Int)\n\nCreate a new pool by breaking up the DNA fragments in an input pool.\n\nThis method breaks up a DNA molecule in a pool p, such that the average length of the fragments is approximately meansize.\n\nIt fragments a molecule by scattering an appropriate number of breakpoints across the molecule, before cutting the molecule at those breakpoints.\n\nnote: Note\nBreakpoints are scattered entirely at random across a molecule. No two or more breakpoints can fall in exactly the same place, as those positions are sampled without replacement.\n\nnote: Note\nThe appropriate number of breakpoints to scatter across a molecule is calculated as:fracLS - 1Where L is the length of the molecule being fragmented, and S is the desired expected fragment size. This calculation assumes breakpoints fall randomly across the molecule (see above note).\n\nnote: Note\nIf a DNA molecule being fragmented is smaller than the desired meansize, then it will not be broken, it will simply be included in the new pool.\n\n\n\n\n\n"
},

{
    "location": "api/pool/#Fragment-1",
    "page": "Molecule Pool",
    "title": "Fragment",
    "category": "section",
    "text": "fragment(p::Pseudoseq.MoleculePool, meansize::Int)"
},

{
    "location": "api/pool/#Pseudoseq.subsample-Tuple{Pseudoseq.MoleculePool,Int64}",
    "page": "Molecule Pool",
    "title": "Pseudoseq.subsample",
    "category": "method",
    "text": "subsample(p::MoleculePool, n::Int)\n\nCreate a new pool by sampling an input pool.\n\nnote: Note\nDNA molecules in the input pool p are selected according to the uniform distribution; no one molecule is more or less likely to be selected than another.\n\nnote: Note\nSampling is done without replacement, so it is impossible for the new pool that is created to recieve one molecule the input pool twice.\n\n\n\n\n\n"
},

{
    "location": "api/pool/#Subsampling-1",
    "page": "Molecule Pool",
    "title": "Subsampling",
    "category": "section",
    "text": "subsample(p::Pseudoseq.MoleculePool, n::Int)"
},

{
    "location": "api/pool/#Pseudoseq.tag-Tuple{Pseudoseq.MoleculePool,Int64}",
    "page": "Molecule Pool",
    "title": "Pseudoseq.tag",
    "category": "method",
    "text": "tag(u::MoleculePool, ntags::Int)\n\nCreate a pool of tagged DNA molecules from some input pool.\n\nThe new tagged pool has the same DNA molecules as the input pool. However, each DNA molecule in the new tagged pool will be assigned a tag in the range of 1:ntags.\n\nFor any tagged molecules in a pool, any other molecules that are derived from that tagged molecule will inherit the same tag. For example, if a DNA fragment in a pool is tagged, and then it is subsequently fragmented during a fragment transform, then all the smaller fragments derived from that long fragment will inherit that long fragment\'s tag.\n\nnote: Note\nWhich fragment gets a certain tag is completely random. It is possible for two distinct DNA molecules in a pool to be assigned the same tag. The likelihood of that happening depends on the size of the tag pool (ntags), and the number of fragments in the pool.\n\n\n\n\n\n"
},

{
    "location": "api/pool/#Tagging-1",
    "page": "Molecule Pool",
    "title": "Tagging",
    "category": "section",
    "text": "tag(p::Pseudoseq.MoleculePool, ntags::Int)"
},

{
    "location": "api/pool/#Flipping-1",
    "page": "Molecule Pool",
    "title": "Flipping",
    "category": "section",
    "text": "flip(p::Pseudoseq.MoleculePool)"
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
    "location": "api/reads/#Exported-functions-1",
    "page": "Reads",
    "title": "Exported functions",
    "category": "section",
    "text": ""
},

{
    "location": "api/reads/#Pseudoseq.make_reads",
    "page": "Reads",
    "title": "Pseudoseq.make_reads",
    "category": "function",
    "text": "make_reads(::Type{PairedEnd}, p::MoleculePool, flen::Int, rlen::Int = flen)\n\nCreate a set of paired-end reads from a pool of DNA molecules p.\n\nflen sets the length of forward read, and rlen sets the length of the reverse read. If you only provide flen, then the function sets rlen = flen.\n\nnote: Note\nIf a molecule in the pool is not long enough to create a forward and/or reverse read, then that molecule will simply be skipped. \n\n\n\n\n\nmake_reads(::Type{SingleEnd}, p::MoleculePool, len::Int)\n\nCreate a set of single-end reads from a pool of DNA molecules p.\n\nlen sets the length of the reads.\n\nThe end (strand) from which the reading begins for each DNA molecule in the pool is determined at random for each molecule, with 50:50 probability.\n\nIf you don\'t provide a value for len, then the function will read each DNA molecule in it\'s entirety.\n\nnote: Note\nIf a molecule in the pool is not long enough to create a forward and/or reverse read, then that molecule will simply be skipped. \n\n\n\n\n\nmake_reads(::Type{TaggedPairs}, p::MoleculePool, flen::Int, rlen::Int = flen)\n\nCreate a set of tagged paired-end reads from a pool of DNA molecules p.\n\nflen sets the length of forward read, and rlen sets the length of the reverse read. If you only provide flen, then the function sets rlen = flen.\n\nWhen a set of TaggedPairs is written to file, the tag information is contained in the R1 read of each read-pair. The first 16bp of each R1 read is a sequence that is the tag, and a following 7bp are a buffer between the 16bp tag, and the rest of the read sequence.\n\nnote: Note\nIf a molecule in the pool is not long enough to create a forward and/or reverse read, then that molecule will simply be skipped. \n\n\n\n\n\n"
},

{
    "location": "api/reads/#Making-reads-1",
    "page": "Reads",
    "title": "Making reads",
    "category": "section",
    "text": "make_reads"
},

{
    "location": "api/reads/#Pseudoseq.mark_errors",
    "page": "Reads",
    "title": "Pseudoseq.mark_errors",
    "category": "function",
    "text": "mark_errors(reads::Reads, rate::Float64)\n\nCreate a new set of reads, with errors, from an input set of reads.\n\nWhen you first create a set of reads using a make_reads method, all the reads in that set are perfect reads. In real sequencing experiments, sequencers make errors when reading a DNA molecule, and are characterised by an error rate. This function lets you simulate this characteristic of sequencers, by marking (at random) positions in the reads that are destined to be errors in the output FASTQ.\n\nnote: Note\nCurrently, every position in every read is equally likely to be marked as an error.\n\nnote: Note\nThe number of errors introduced into the reads can be calculated.E = NRWhere E is the number of errors, N is the total number of bases in your set of reads, and R is the error rate you provide to this function.\n\n\n\n\n\n"
},

{
    "location": "api/reads/#Introducing-errors-1",
    "page": "Reads",
    "title": "Introducing errors",
    "category": "section",
    "text": "mark_errors"
},

{
    "location": "api/reads/#Pseudoseq.generate",
    "page": "Reads",
    "title": "Pseudoseq.generate",
    "category": "function",
    "text": "generate(filename::String, reads::Reads)\n\nWrite the reads out to a FASTQ formatted file with the given filename.\n\nIf this method is used with a paired-end read type, then the FASTQ file will be interleaved; all R1 reads will be odd records, and all R2 reads will be even records in the file.\n\nnote: Note\nReads are named according to the sequence in the input genome they came from. e.g. @Reference_1_R1 means the first sequence in the genome, and @Reference_2_R1 means the second sequence in the genome.\n\n\n\n\n\ngenerate(R1name::String, R2name::String, reads::Reads{<:PairedReads})\n\nThis method only works for paired reads. Instead of interleaving R1 and R2 reads in a single FASTQ file, R1 and R2 reads are partitioned into two seperate FASTQ files.\n\n\n\n\n\n"
},

{
    "location": "api/reads/#Generating-FASTQ-files-1",
    "page": "Reads",
    "title": "Generating FASTQ files",
    "category": "section",
    "text": "generate"
},

{
    "location": "api/chromosome-blueprint/#",
    "page": "Build-a-Genome",
    "title": "Build-a-Genome",
    "category": "page",
    "text": ""
},

{
    "location": "api/chromosome-blueprint/#API:-Build-a-Genome-1",
    "page": "Build-a-Genome",
    "title": "API: Build-a-Genome",
    "category": "section",
    "text": ""
},

{
    "location": "api/chromosome-blueprint/#Exported-functions-1",
    "page": "Build-a-Genome",
    "title": "Exported functions",
    "category": "section",
    "text": ""
},

{
    "location": "api/chromosome-blueprint/#Pseudoseq.plan_chrom",
    "page": "Build-a-Genome",
    "title": "Pseudoseq.plan_chrom",
    "category": "function",
    "text": "Create an empty blueprint for n copies of a chromosome of len base pairs.\n\n\n\n\n\n"
},

{
    "location": "api/chromosome-blueprint/#Making-an-empty-chromosome-blueprint-1",
    "page": "Build-a-Genome",
    "title": "Making an empty chromosome blueprint",
    "category": "section",
    "text": "plan_chrom"
},

{
    "location": "api/chromosome-blueprint/#Pseudoseq.plan_repetition",
    "page": "Build-a-Genome",
    "title": "Pseudoseq.plan_repetition",
    "category": "function",
    "text": "plan_repetition(cb::ChromosomeBlueprint, from::UnitRange{Int}, to::UnitRange{Int})\n\nPlan a repetition in a chromosome, where the bases in the from region of the chromosome, are guaranteed to occur again in the to region of the chromosome.\n\nEvery copy of the chromosome will have this same repetition.\n\nCreates a new chromosome blueprint, based on the input blueprint cb.\n\nnote: Note\nIn the new blueprint, the to region of the planned chromosome will be consumed, and cannot be used to plan any other subsequently added features.\n\n\n\n\n\nplan_repetition(cb::ChromosomeBlueprint, from::Int, to::Int, size::Int)\n\nPlan a repetition in a chromosome, where the bases in the from:(from + size - 1) region of the chromosome, are guaranteed to occur again in the to:(to + size - 1) region of the chromosome.\n\nEvery copy of the chromosome will have this same repetition.\n\nCreates a new chromosome blueprint, based on the input blueprint cb.\n\n\n\n\n\nplan_repetition(cb::ChromosomeBlueprint, intervals::Vector{UnitRange{Int}})\n\nA conveinience method of plan_repetition. Designed to ease the process of planing a series of repetitions in a chromosome.\n\nEvery copy of the chromosome will have these same repetitions.\n\nCreates a new chromosome blueprint, based on the input blueprint cb.\n\ntip: Tip\nUse the suggest_regions function to help decide on a set of sites to make heterozygous.\n\nnote: Note\nThe number of intervals provided must be an even number. This is because intervals 1 & 2  define the first repeat, intervals 3 & 4 define the second, and so on.\n\nnote: Note\nIn the new blueprint, the regions the repetitions occupy, have been consumed, and cannot be used to plan any other subsequently added features.\n\n\n\n\n\n"
},

{
    "location": "api/chromosome-blueprint/#Planning-motif-repetitions-along-chromosomes-1",
    "page": "Build-a-Genome",
    "title": "Planning motif repetitions along chromosomes",
    "category": "section",
    "text": "plan_repetition"
},

{
    "location": "api/chromosome-blueprint/#Pseudoseq.plan_het",
    "page": "Build-a-Genome",
    "title": "Pseudoseq.plan_het",
    "category": "function",
    "text": "plan_het(cb::ChromosomeBlueprint, pos::Int, alleles::Vector{DNA})\n\nPlan heterozygosity between copies of a chromosome at position pos.\n\nThe alleles vector must contain a nucleotide for each copy of the chromosome.\n\nFor example if you provided a vector of [DNA_A, DNA_C], for the heterozygous site you\'re defining at pos, the first copy of the chromosome will have an A, and the second copy of the chromosome will have a C.\n\nCreates a new chromosome blueprint, based on the input blueprint cb.\n\ntip: Tip\nUse the suggest_alleles function to help decide on a set of alleles to use.\n\nnote: Note\nIn the new blueprint, the positions that were used to plan the heterozygosity, have been consumed, and cannot be used to plan any other subsequently added features.\n\n\n\n\n\nplan_het(cb::ChromosomeBlueprint, pos, alle...)\n\nA generic method of plan_het.\n\nCreates a new chromosome blueprint, based on the input blueprint cb.\n\nThe pos argument determines which sites in a chromosome are made heterozygous, and depending on the type of argument provided, the behaviour differs slightly:\n\nIf pos is an integer, then that many available sites are selected at random to be heterozygous.\nIf pos is a float, then it is treated as a proportion of the length of the, chromosome, which is converted into an integer, and that many sites are selected at random to be heterozygous.\nIf pos is a vector of integers, it is trated as a list of positions the user wants to be heterozygous.\n\nThe alle... argument is a vararg argument which defines the how bases will be  allocated at each heterozygous site.\n\nDepending on the type of the argument provided, the behaviour differs slightly, based on the different methods of suggest_alleles:\n\nIf alle... is a single integer value. It is taken to mean the number of states each heterozygous site has. At a minumum this number must be two, as a site with only one state is not heterozygous by definition. For each site, which chromosome copies get which state is determined at random.\nIf alle... is a vector of integers (one integer for every chromosome copy), then the nth value in the vector dictates the group the nth chromosome copy belongs to. For example for a triploid, the vector [2, 1, 2] means that there are two groups (group 1 and group 2), and that the second chromosome copy of the three belongs to group 1, and the other two copies belong to group\nChromosome copies that belong to the same group will get the same base at\neach heterozygous site. In our example, copies 1 and 3 will get the same base, as they have been given the same group number. Chromosome copy 2 will end up with a different base. Which base corresponds to which group number will be determined at random.\nIf alle... is multiple vectors of integers, then each vector constitutes a group, and each value in the vector determines which copies are in that group. It is alternative way to giving the same information as you can with a single vector (as above). E.g. consider the example of [2, 1, 2] from the previous point. This can be expressed as two vectors of [2] and [1, 3]: The first vector contains a number indicating the second chromosome copy, and the second vector contains numbers representing the first and third chromosome copy.\nIf alle... is a vector of nucleotide vectors, then the nth vector of nucleotides determines which chromosome copy recieves which base at the nth heterozygous position. For example, for a triploid, [[DNA_A, DNA_A, DNA_T], [DNA_G, DNA_C, DNA_C]] would mean that at the first heterozygous site, the first and second copies of the chromosome would recieve an DNA_A base, and the third copy would recieve a DNA_T base. For the second heterozygous position, the first copy of the chromosome would recieve a DNAG base, and the other two copies will recieve a `DNAC` base.\n\n(See also: plan_repetition)\n\n\n\n\n\n"
},

{
    "location": "api/chromosome-blueprint/#Planning-heterozygosity-between-chromosome-copies-1",
    "page": "Build-a-Genome",
    "title": "Planning heterozygosity between chromosome copies",
    "category": "section",
    "text": "plan_het"
},

{
    "location": "api/chromosome-blueprint/#Pseudoseq.fabricate",
    "page": "Build-a-Genome",
    "title": "Pseudoseq.fabricate",
    "category": "function",
    "text": "fabricate(file::String, cb...)\n\nFabricate the sequence(s) planned in a number chromosome blueprints.\n\nThe fabricated sequences will be written out to the FASTA formatted file file.\n\ncb... should be provided as a series of blueprint, seed-sequence pairs. E.g. fabricate(\"mygenome.fasta\", chr1plan => ch1seq, ch2plan => ch2seq).\n\n\n\n\n\nfabricate(fw::FASTA.Writer, cb...)\n\nFabricate the sequence(s) planned in a number chromosome blueprints.\n\nThe fabricated sequences will be written out to the FASTA formatted file fw.\n\ncb... should be provided as a series of blueprint, seed-sequence pairs. E.g. fabricate(\"mygenome.fasta\", chr1plan => ch1seq, ch2plan => ch2seq).\n\n\n\n\n\nfabricate(cb::ChromosomeBlueprint)\n\nFabricate the sequence(s) planned in the chromosome blueprint.\n\nA random DNA sequence will be generated to use as a seed sequence.\n\nA sequence will be built for each chromosome copy in the blueprint.\n\n\n\n\n\nfabricate(cb::ChromosomeBlueprint, seed::BioSequence{DNAAlphabet{2}})\n\nFabricate a DNA sequence by applying the planned features in a chromosome blueprint, to some initial starting seed sequence.\n\nA sequence will be built for each chromosome copy in the blueprint.\n\n\n\n\n\n"
},

{
    "location": "api/chromosome-blueprint/#Fabricate-1",
    "page": "Build-a-Genome",
    "title": "Fabricate",
    "category": "section",
    "text": "fabricate"
},

{
    "location": "api/chromosome-blueprint/#Pseudoseq.suggest_regions",
    "page": "Build-a-Genome",
    "title": "Pseudoseq.suggest_regions",
    "category": "function",
    "text": "suggest_regions(cp::ChromosomeBlueprint, size::Int, n::Int)\n\nA useful utility function to assist planning features in a chromosome blueprint.\n\nThis method returns a vector of non-overlapping, regions of the chromosome planned represented by the ChromosomeBlueprint cb.\n\nThese regions are free regions: they are untouched by any other planned features, and so may be used when planning other features (See also: plan_repetition, plan_het).\n\nThe regions will be sizebp in length.\n\nwarning: Warning\nThis function was designed for use interactively in a julia session.If this method cannot find n free regions of the size you\'ve asked for, it will still give an output vector containing the regions it did manage to find, but it will issue a warning to the terminal. This will get increasingly likely as you fill the chromosome blueprint up with features.If you are using this method in a script or program where you depend on a reliable number of regions, either add a check to make sure you got the number of regions you need, or use this method interactively, and hard code an appropriate output into your script.\n\n\n\n\n\n"
},

{
    "location": "api/chromosome-blueprint/#Pseudoseq.suggest_alleles",
    "page": "Build-a-Genome",
    "title": "Pseudoseq.suggest_alleles",
    "category": "function",
    "text": "suggest_alleles(groups::Vector{Int})\n\nSuggest an allele pattern for a single heterozygous site.\n\nThe groups vector should have a number of  The nth value in the groups designates a group to the nth chromosome copy.\n\nE.g. consider for a triploid, the vector [2, 1, 2] means that there are two groups (group 1 and group 2), and that the second chromosome copy of the three belongs to group 1, and the other two copies belong to group 2.\n\nChromosome copies that belong to the same group will get the same base at each heterozygous site.\n\nFor the vector [2, 1, 2], this method might return a suggested allele pattern such as [DNA_A, DNA_G, DNA_A].\n\nnote: Note\nWhich base corresponds to which group number will be determined at random.\n\n\n\n\n\nsuggest_alleles(groups::Vector{Int}...)\n\nSuggest an allele pattern for a single heterozygous site.\n\nMultiple integer vectors make up groups.\n\nEach vector defines a group, and each value in a vector determines which copies are allocated to that group.\n\nE.g. consider for a triploid, two vectors of [2] and [1, 3] passed as arguments to groups. They define two groups, the first vector contains a number denoting the second chromosome copy, and he second vector contains numbers representing the first and third chromosome copy.\n\nFor the vectors [2] and [1, 3], this method might return a suggested allele pattern such as [DNA_A, DNA_G, DNA_A].\n\nnote: Note\nWhich base corresponds to which group number will be determined at random.\n\n\n\n\n\nsuggest_alleles(ncopies::Int, ngroups::Int)\n\nSuggest an allele pattern for a single heterozygous site.\n\nBy providing a number representing the number of copies of a chromosome, and a number representing the number of groups (which at a minimum must be 2). This method of suggest_alleles will randomly allocate the ncopies of the chromosome into the ngroups groups.\n\nFor example, if you ask for a suggested allele pattern for 3 copies of a chromosome, and 2 groups, this method might return a suggested pattern of [DNA_G, DNA_A, DNA_G].\n\nnote: Note\nWhich base corresponds to which group number will be determined at random.\n\n\n\n\n\nsuggest_alleles(npositions::Int, args...)\n\nSuggest allele patterns for numerous heterozygous sites.\n\nThis method of suggest_alleles calls another suggest_alleles with args... npositions times, and returns a vector of the results of each call.\n\n\n\n\n\n"
},

{
    "location": "api/chromosome-blueprint/#Utility-functions-1",
    "page": "Build-a-Genome",
    "title": "Utility functions",
    "category": "section",
    "text": "suggest_regions\nsuggest_alleles"
},

]}
