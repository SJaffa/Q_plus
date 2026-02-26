# Record of a specific analysis [[raw](https://gitlab.com/rsspdc/checklists/-/raw/master/checklists/rsspdc-specific-record.md)]

version: alpha

*[Research Software Sharing, Publication, & Distribution Checklist](https://gitlab.com/rsspdc/checklists)*

Considerations for publishing code which runs a specific analysis that underpins some result to be published in the academic literature.
The emphasis here is on making the work narrowly reproducible i.e. the analysis of the same data can produce the same result when it is re-run.
This is a starting point for making results robust (different analysis, same data) and replicable (same analysis, different data), and ultimately generalisable (different analysis, different data).
The other emphasis is on making the work 'verifiable', exposing the complete step-wise detail of the reasoning underpinning the analysis so that it can be scrutinised and understood.

## 📒Source control

*How can you keep track of the history of your project and collaborate on it?*

- [ ] Uses git (or other source control tool)
	- [ ]  	🥉Bronze *(Easy)*: Using version control but has a shallow project history, just placed in git for distribution
	- [ ]  	🥈Silver *(Intermediate)*: Longer project history, commit messages of mixed quality, some large messy changes
	- [ ]  	🥇Gold *(Hard)*: Silver plus - Well written commit messages, nice granular commits making discrete self-contained changes. Tags, releases, or branches at major project milestones, maybe some contributions from other users
	- [ ]  	🏆Platinum *(MAXIMUM OVERKILL)*: Gold plus - Some from: [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/); Clean history with a consistent rebasing/merging strategy; Signed commits from all contributors; Contributions go through a consistent workflow like, issues, then a pull request from a branch.

<details>
<ul>
	<li>
	Whilst you can simply use git and a git forge / git hosting service as a way of distributing your project doing so misses out on a lot of benefits of using git as a part of your workflow from the beginning of your project.
	A well maintained git history is much like a well kept lab notebook for a data analysis project.
	Well authored commit messages detail why you changed what you changed, providing context for the development of the project.
	You can by 'checking out' a commit open a window onto any point in the history of your project, of which you took a snapshot by making a commit.
	You can collaborate on your project with other git users asynchronously, it can be a great tool for distributed collaborative authorship not just of software but also of prose.
	A good example of this is [The Turing Way: a how to guide for reproducible data science](https://github.com/the-turing-way/the-turing-way/).
	A nice read making this case is [Not just for programmers: How GitHub can accelerate collaborative and reproducible research in ecology and evolution](https://doi.org/11.1111/2041-210X.14108), though I would carefully weigh the long term issues often created by building on proprietary infrastructure like github as opposed to using an open git forge instead.
	</li>
	<li>
	Interaction with data and other big non-text / binary files.
	</li>
	<ul>
		<li>
		git is not good at managing changes to large binary files like images, all versions of the file will be retained in the git history which can make it grow fast if they change frequently.
		</li>
		<li>
		You can exclude large datafiles from tracking by git by adding them, or a pattern matching them such as: `data/*` (ignore all files in in the data folder) or `*.data` (ignore all files with that end in `.data`) to the `[.gitignore](https://git-scm.com/docs/gitignore)` file.
		</li>
		<li>
		There are a variety of different tools which can be used to manage versioning large binary objects in a git like fashion and which integrate with a git based workflow.
		If you use one, and which you choose, may depend on your specific needs.
		Examples: [git-annex](https://git-annex.branchable.com/), [git LFS](https://git-lfs.com/), [Data Version Control (DVC)](https://dvc.org/), [lakeFS](https://lakefs.io/)
		This is generally most relevant for intermediate data objects which are of potential interest to downstream users of the processed data objects.
		for example machine learning model weights from different training runs, annotated single cell sequence count matrices in things like seurat objects.
		If you have deposited your raw data into a public repository it does not need to be duplicated in such a system indefinietly but it might be useful whilst working on the project to have it 'cached' in one.
		If possible start by importing the data & metadata into your own project from it's public repository as a test of FAIRness of the data.
		This also means that anyone using your dataset has a clear example of how to import it into a working environment.
		</li>
	<li>
	Signed and Timestamped git commits
	</li>
	<ul>
		<li> 
		It is possible to cryptographically sign your git commits, this can be used to increase confidence that you are the author of a signed commit as someone would have to compromise your private key in order to impersonate you.
		If your key is part of a [web of trust](https://en.wikipedia.org/wiki/Web_of_trust) or other (public key infrastructure (PKI)](https://en.wikipedia.org/wiki/Public_key_infrastructure) people can see that other people / institutions attest that person with this key is who they say they are.
		Whilst typically used for things like helping to protect critical open source infrastructure from supply chain attacks signed commits on academic code bases could be used to provide additional provenance information.
		</li>
		<li> 
		It is sometimes desirable for git commits to demonstrably have been signed at a given time, your system time is recorded by default but this can be trivially spoofed.
		The [opentimestamps](https://opentimestamps.org/) protocol which can be used to generate cryptographic attestations to the time at which a commit was made.
		This aims to establish a lower bound on how long ago the committed code was authored.
		The [opentimestamps-client](https://github.com/opentimestamps/opentimestamps-client) [integrates with git](https://github.com/opentimestamps/opentimestamps-client/blob/master/doc/git-integration.md) to provide timestamps for individual [GPG](https://gnupg.org/) signed commits
		</li>
	</ul>
</ul>
</details>

## ©Licencing

*On what terms can others use your code, and how can you communicate this?*

- [ ] Project is suitably licensed
	- [ ]  	🥉Bronze *(easy)*: There is a LICENSE file in the repository for a license which meets one of the [OSI](https://opensource.org/osd), [Debian](https://www.debian.org/intro/free), or [FSF/GNU](https://www.gnu.org/philosophy/free-sw.en.html) definitions of free/libre or open source software.
	Or for any contents that are not software a [Creative Commons](https://creativecommons.org/) or similar license which meets the [free cultural works](https://freedomdefined.org/Licenses) definition.
	- [ ]  	🥈Silver *(easy)*: If any prose/documentation or images are licenced differently from the code in the project this is indicated and those licences provided.
	If licences have an attribution requirement there is easy to copy text/links for appropriate attribution.
	- [ ]  	🥇Gold *(intermediate)*: Uses [Software Package Data Exchange (SPDX)](https://spdx.org/licenses/) license identifiers for every file/suitable unit of code.
	With a tool such as [REUSE.software](https://reuse.software/) to automate and standardise the process.
	- [ ]  	🏆Platinum *(intermediate)*: all previous tiers plus any images have licensing information embedded in their metadata.

<details>
<ul>
	<li>
	Once you have selected you license include a plain text copy of it in the root of your repository in a file named LICENSE.
	Plain text versions of popular open license are widely available.
	Files with this name are often identified by software forges and a link to them created on the repository home page when they are present.
	You may need to include your name and the date in your copy of the license file if indicated.
	</li>
	<li>
	The content of a repo of this form is generally a mixture of code, images (often graphs), data, and prose.
	In this context it may be preferable to have licenses for code, prose, and other assets such as graphs. e.g. all code under a GPLv3 license, and all images, prose and datasets, under a CC BY-SA license.
	</li>
	<li>
	**All software needs a license if you want to permit others to reuse it.**
	It is important to give some thought to the type of license which best suits your project, it is a choice which can have significant long term implications.
	Checkout [the turing way chapter on licensing for an introduction to the subject](https://the-turing-way.netlify.app/reproducible-research/licensing.html).
	If you have no time some pretty safe choices are: For a permissive license, the Apache 2.0.
	This would allow the re-use of your work in closed commercial code.
	For a 'copyleft' license, the GPLv3 (AGPL for server-side apps).
	This requires that anyone distributing software containing your code or derivatives of it share the source code with the people they distributed it to.
	</li>
	<li>
	If you are including external code in your package then you should check that their licenses are compatible and you are legally allowed to distribute your code together in this way.
	Checkout this [resource on license compatibility](https://the-turing-way.netlify.app/reproducible-research/licensing/licensing-compatibility.html).
	</li>
	<li>
	If you have a big project with a lot on differently licensed content or need a standard way to provide license information for binary assets you might want to checkout the [REUSE.software](https://reuse.software/) tool, you can also embed licensing information directly into image metadata to help ensure that it remains associated with the image when/if it is reused, [IPTC](https://www.iptc.org/std/photometadata/specification/IPTC-PhotoMetadata#copyright-notice) provides an extended image metadata schema with a copyright notice field which could house an [SPDX](https://spdx.dev/) licence code.
	It also allow for the direct embedding of image alt text.
	</li>
</ul>
</details>

## 📖Documentation

*How do people know what your project is, how to use it and how to contribute?*

- [ ] Project has suitable documentation 
	- [ ] 	🥉Bronze *(easy)*: Project has README file that either contains directly or links to resources which answer these questions:
		- [ ] Provides a description of the project structure so that the user knows which directories to find things in, possibly including a visual representation of the structure
		- [ ] Contains instructions with sufficient detail for someone else to re-run the analysis
		- [ ] What is the name or title of the project for which this code is record?
		- [ ] Why did you peform this analysis?
		- [ ] What problem does it address and how?
		- [ ] How do I set up an environment to re-run it? (with example)
		- [ ] What inputs do I need?
		- [ ] What outputs should I expect?
		- [ ] How do I cite the project?
		- [ ] Who contributed to the project and what did they contribute?
		- [ ] Who should I contact, about what, and how?
		- [ ] What should I expect / Not expect if I contact you?
		- [ ] How is the project Licensed?
		- [ ] Table of Contents (if long enough to benefit from one)
		- [ ] Bibliography (if reference is made to external resources)
	- [ ]  	🥈Silver *(intermediate)*: Documentation / comments explain WHY things are done in the code
	- [ ]  	🥇Gold *(intermediate)*: You have a simple worked example of your analysis methods with example data to illustrate its soundness in a simple case
	- [ ]  	🏆Platinum *(MAXIMUM OVERKILL)*: The authors manuscript of you paper is a literate programming artefact with numbers, tables and figures programaticaly generated and built in a reproducible computational environment.
			Additional documentation of details not suitable for the main manuscript are in the supplementary material also as literate programming artefacts, more lengthy and complex elements of the analysis might be run as a pipeline or scripts then read from for inclusion in the manuscript.
			Building the manuscript (without caching) reruns the complete analysis and produces the same output (you may need to set random seeds for any analysis that makes use of pseudorandom number generation).

In some cases for small and simple projects a README file is sufficient documentation.
This may genuinely be all that you need, or inadequate to the task depending on the project.

<details>
<ul>
	<li>The rationale is key for code pertaining to a specific analysis so that the intent and reasoning of the author is clear to the reader</li>
	<li>
	literate programming notebook tools like Jupyter and Rmarkdown/Quarto are a great tool for this sort of code output as they permit you to contextualise the choices made during an analysis, visualise, and interpret your results all in the same place.
	They can also provide robust provenance of results.
	It can be made clear that this data, was analysed by this code, in this compute environment, and produced these outputs as the output is built by executing the notebook.
	</li>
	<li>
	Analyses may be sufficiently complex that you don't want all of them in a literate programming form especially if you are using literate programming tools to author a conventional manuscript.
	You can treat the full details of your analysis as a supplementary method and use a combination of more conventional software documentation tools and literate programming to provide a full detail version of your analysis then reference these objects in your manuscript.
	If working in R the {[targets](https://books.ropensci.org/targets/)} tool can be a nice way of managing a workflow like this as you can cache results like graphs that might be complex and expensive to generate then reference them succinctly in your manuscript document.
	</li>
	<li>
	A good test of reproducibility of your analysis is to use CI/CD tools to build your outputs from source code and never commit the outputs themselves to the repo, you can then serve these build artefacts as a static website of your manuscript and any other documentation from your project.
	You can also bundle them with the rest of the code when creating a snapshot of the project to archive on a platform like [zenodo](https://zenodo.org/).
	Obviously it is best to avoid re-running lengthy computationally intensive analyses every time you push to a repository.
	You can avoid this by caching results and/or only triggering rebuilds when a commit is tagged a certain way, for example a new version number.
	This way it may still be possible to run your more computationally intensive code on your build system, if your build system has adequate computational resources.
	</li>
	<li>
	If you release multiple versions of a particular analysis it can be a good idea to include a CHANGELOG file in your project documenting things which have changed since the previous version.
	</li>
</ul>
</details>

## 🔗Making Citable

*How should people make reference to your project and credit your work?*

- [ ] Record is Citable
	- [ ]  	🥉Bronze *(easy)*: A [CITATION.cff](https://citation-file-format.github.io/) file exists in the code repository to provide citational metadata about your project
	- [ ]  	🥈Silver *(easy)*: bronze plus the project has persistent resolvable identifier such as a DOI or SWHID, with which it can be referenced, which has been minted for the project using a tool like [zenodo](https://zenodo.org/) or [Software Heritage's Archive](https://www.softwareheritage.org/) to store an archival copy of the project.
	- [ ]  	🥇Gold *(intermediate)*: silver plus:
		- [ ] Contributions are credited using a suitable contributor roles ontology or taxonomy (CROT) such as [CrediT](https://credit.niso.org/), [ScoRo](http://www.sparontologies.net/ontologies/scoro), [CRO](https://github.com/data2health/contributor-role-ontology), or [TaDiRAH](https://tadirah.info/).
		- [ ] All contributors are identified by their [ORCID](https://orcid.org/) or other suitable persistent identifier
	- [ ]  	🏆Platinum *(mixed)*: gold plus any two or more from:
		- [ ] (easy) All research institutions are identified by their [ROR](https://ror.org/) ID
		- [ ] (intermediate) Versioned persistent identifier with automation to update snapshots on zenodo or similar tool when a new version is created.
		- [ ] (intermediate) Annotating work cited in this work with the [Citation Typing Ontology (CiTO)](http://purl.org/spar/cito)
		- [ ] (hard) Your environment is defined with Nix or Guix - this might not seem like it contributes to making software more citable see details below for why this is the case.

<details>
<ul>
	<li>
	**The code underpinning a publication is a part of your methods.**
	The legacy publishing system lacks a suitable and convenient way of including your analysis code as a part of your methods section so to work around this it is best to make your code a seperate citable object and just cite it in your methods section.
	</li>
	<li>
	Including a [CITATION.cff](https://citation-file-format.github.io/)(Citation File Format) file in your project repo is a simple way of making your code citable.
	The format is readable in YAML and permits the provision of the metadata needed for citation.
	</li>
	<li>
	[Zenodo](https://zenodo.org/) permits you to mint a [digital object identifier (DOI)](https://www.doi.org/) for your code and makes a snapshot of it, importing citational metadata from a [CITATION.cff](https://citation-file-format.github.io/) file or a .zenodo.json file.
	This makes it persistently identifiable and easy to integrate with citation management tools which can import the citation metadata given a DOI.
	</li>
	<li>
		Nix and GUIX
		<ul>
			<li>
				General software repositories may not make specific provision for citation of software packages in the academic fashion.
				However some provide, what is for some use cases, a superior form of 'citation' of their own sources i.e. a complete 'software bill of materials (SBOM)'.
				This is a list of all the code used in another piece of code, its dependencies, and their dependencies recursively, along with all of their versions.
				For example [Nix](https://nixos.org/) can do this but [Guix](https://guix.gnu.org/) has perhaps the most comprehensive in its approach.
				It not only provides all information necessary for a complete 'SBOM' but, it can [bootstrap](https://guix.gnu.org/en/manual/en/html_node/Bootstrapping.html) software packages in its repository from source with an extremely minimal fixed set of binaries, an important capability for [creating somewhat trustworthy builds](https://www.cs.cmu.edu/~rdriley/487/papers/Thompson_1984_ReflectionsonTrustingTrust.pdf).
				This creates a compute environment which is not only reproducible but verifiable, meaning the source of all of an environment's dependencies can in theory be scrutinised.
				It also adopts an approach to commit signing and authorisation of signers that gives it a [currently uniquely complete supply chain security architecture](https://doi.org/10.22152/programming-journal.org/2023/7/1).
				Packages or 'derivations' are 'pure functions' in the sense that only their inputs effect their outputs and they have no side-effects, package builds are sandboxed to prevent dependencies on any external source not explicitly provided as an input and inputs are hashed to ensure that they cannot differ for the value expected when they were packaged.
				This gives these technologies an unrivaled ability to readily demonstrate the reproducibility and provenance of compute environments specified using them.
			</li>
			<li>
				Whilst not yet full implemented and adopted these technologies also aford some fascinating opertunities for seemless access to archival versions of software in the future.
				Due to the similarities in the content based addressing used by Git, Nix, Guix, IPFS (Interplanetary file System) and software heritage's IDs it may be possible to construct an approach to archiving, distributing and caching sources of packages in a way that would ensure that low demand archived software sources and high demand current packages can be distributed transparently through the same mechanism.
				This would in theory permit the reconstruction of any historically specified compute environment that had been archived with no changes to normal workflow, other than perhaps a longer build time.
				This approach also makes the creation of 'mirrors' of the archive relatively simple and requires no client side changes as an IPFS resource will be resolved irrespective of the node on which it is stored.
				See: [NLnet Software heritage and IPFS](https://nlnet.nl/project/SoftwareHeritage-P2P/), [Tweag - software heritage and Nixpkgs](https://www.tweag.io/blog/2020-06-18-software-heritage/), [John Ericson - Nix x IPFS Gets a New Friend: SWH (SoN2022 - public lecture series)](https://www.youtube.com/watch?v=DjJyPzwEzmU)
			</li>
		</ul>
	</li>
</ul>
</details>

## ✅Testing

*How can you test your project so you can be confident it does what you think it does?*

- [ ] Project has undergone suitable testing
	- [ ]  	🥉Bronze *(easy)*: Includes a minimal test data set necessary to demonstrate the basic functionality of the analysis
	- [ ]  	🥈Silver *(easy)*: Includes test datasets which cover a range of outcomes of the analysis
	- [ ]  	🥇Gold *(intermediate)*: You are using unit tests and an automated testing framework to check the correctness of core steps of your analysis
	- [ ]  	🏆Platinum *(hard)*: Have your analysis code written and tested on preliminary or simulated data in advance of recieving your principle dataset with a copy of your code from this time archived and referenced in a pre-registration or registered report.

<details>
<ul>
	<li>
	Whilst you can make use of unit tests / automated testing frameworks in this context (see the software packages checklist testing section for more) it is not always the best fit.
	A very good thing to do is to have an example dataset that you can perform your analysis on that is different from your new data.
	If you are testing a hypothesis it's nice to have test datasets which simulate rejecting and accepting your null hypothesis.
	</li>
	<li>
	If you write, test, have reviewed, and publically deposit, the code that you plan to run to test the primary/pre-planned outcomes in your study and are able to run this same code on your data unchanged this inspires confidence in your testing regimend.
	Code associcated with incedental finding which may be suggestive for hypotheses to test in future work obviously cannot be subject to pre-registration but other best practices can be followed.
	</li>
	<li>
	When using real data or downsampled real data be sure that your testing covers any edge cases that may not have arrising enyour example data.
	</li>
	<li>
	When using simulated data be sure to include the method by which you simulated the data and any random seeds which may be needed to re-generate it.
	</li>
</ul>
</details>

## 🤖 **Automation**

*What tasks can you automate to increase consistency and reduce manual work?*

- [ ] Suitable automations are in place
	- [ ]  	🥉Bronze *(easy)*: 1 from this list of processes are automated
		- use of an environment management tool
		- use of a literate programming / computational notebook
		- use of a pipeline manger or make-like tool
		- use of a linter / formatter
		- use of continious integration / continious deployment
		- use of git hooks
		- automated minting of new persistent identifiers on release tagging
		- ...
	- [ ]  	🥈Silver *(easy)*: 2-3 from the above list of processes are automated
	- [ ]  	🥇Gold *(intermediate)*: 4+ from the above list of processes are automated
	- [ ]  	🏆Platinum *(hard)*: note that difficuly is somewhat project dependent
		Your manuscript and its supplements are generated and served on a website after being built from your CI/CD pipeline.
		All statistics and data visualisations in your manuscript are generated programatically by your analysis pipeline from your raw data in CI/CD.
		Results are cached such that if you, for example, change the formatting of a graph only the plotting and rendering code needs to be re-run, but if you change the data the entire pipeline is re-rerun.

<details>
	<ul>
		<li>

		</li>
	</ul>
</details>

## 👥Peer review / Code Review

*How can you get third party endorsement of and expert feedback on your project?*

- [ ] Code has been subject to a review indicating that someone else could re-run the analysis
	- [ ]  	🥉Bronze *(easy)*: Someone other than you has checked over your project, given you feedback and told you they are reasonably confident they could re-run your analysis without your help.
	- [ ]  	🥈Silver *(easy)*: Someone other than you has scessfully re-run your analysis using only your documentation, (preferably in a different compute environment, such as a different computer/compute cluster)
	- [ ]  	🥇Gold *(intermediate)*: You have a review from [CODECHECK](https://codecheck.org.uk/), [ReproHack](https://www.reprohack.org/) or equivalent and have incorporated suggestions for improving reproducibility from these reviews.
	- [ ]  	🏆Platinum *(intermediate)*: You have reviews which go beyond checking the ability to re-run your code but which also review it's technical correctness

<details>
<ul>
	<li>
	Most of the places that offer code peer review are focused on software packages not code that is specific to your analysis.
	This makes sense as reviewer time is fairly scarce so focusing it on code that others are more likely to reuse is reasonable.
	</li>
	<li>
	If your code underpins a publication then in theory it may get reviewed as a part of the regular peer review process although in practice this does not appear to be all that common.
	If the journal to which you are submitting your work has no policy on the code review, and your reviewers do not take an interest in reviewing your code - even just the checking if it runs for them then you may wish to take responibility for the review of this work into your own hands.
	[CODECHECK](https://codecheck.org.uk/) will independently verify that they can run your code, but its correctness is not in their scope.
	</li>
<ul>
</details>

## 📦Distribution

*How can people install or access the software emerging from your project?*

- [ ] Project is distributed in a suitable fashion
	- [ ]  	🥉Bronze *(Easy)*: Code and data (barring privacy related access restrictions) are in public repositories.
	- [ ]  	🥈Silver *(Intermediate)*: Detailed instructions on how to fetch, install and configure the tools and data needed re-run your analysis, and how to re-run the analysis in the described environment.
	- [ ]  	🥇Gold *(Intermediate)*: Project is in a reproducible interactive environment such as those offered by [binder](https://binderhub.readthedocs.io/en/latest/index.html) or [renku](https://renkulab.io/).
	- [ ]  	🏆Platinum *(Hard)*: Gold plus - Your project is built and served as a website using continuous integration and deployment tools such that your analysis is run on your data in a reproducible compute environment and computational results like graphs and statistics are programatically inserted into your output. (It is best to have some form of caching when doing this).

<details>
<ul>
	<li>
	What it means to distribute one off analysis code is somewhat different from distributing a package or pipeline as this goal is different.
	The aims are to share how you did what you did in a research output, and to provide a record of the Provence of your results.
	</li>
	<li>
	Distribution in this context is closer to documentation and literate programming tools like [Jupyter](https://jupyter.org/) Notebooks, [jupyter book](https://jupyterbook.org/en/stable/intro.html) and [Rmarkdown](https://rmarkdown.rstudio.com/)/[Quarto](https://quarto.org/) lend themselves to this task very well.
	You can serve a static website based on the notebooks which perform, explain and interpret your analysis using these tools.
	You can even demonstrate their computational reproduciblity by building these outputs using continuous integration and deployment tools (CI/CD) tools available on gitforges like gitlab and github.
	</li>
	<li>
	Tools like [binder](https://binderhub.readthedocs.io/en/latest/index.html), [renku](https://renkulab.io/) allow you to share your analysis environment so that people can pick up your analysis where you left off in an interactive environment so that they can tweak your code and explore your data as they wish.
	</li>
	<li>
	An excellent way to share the record of a specific analysis is to use all these tools in conjunction.
		<ul>
		<li>
		Perform your analysis in reproducible compute environment specified using a tool like [renku](https://renkulab.io), [binder](https://binderhub.readthedocs.io/en/latest/index.html) or a [Nix flake](https://wiki.nixos.org/wiki/Flakes) which will allow this environment to readily be shared with others.
	Write your manuscript using the literate programming tools and server this as static web page as a way of pre-printing your manuscript.
	Revise it with your collaborators using the issues and pull/merge request features of a git forge.
	If you make repo citable by adding the appropriate metadata and using [zenodo](https://zenodo.org/) to mint a DOI it is as citable as if it were deposited in a pre-print server, but probably looks a lot better.
	Then your entire project history is available in your git history.
	Your computational reproducibility is evidenced by the ability to build the output in the computational environment that you specified to serve your web page with the manuscript.
		</li>
		</ul>
	</li>
	<li>
	Where you have very large and computationally intensive upstream analyses, as is common for example in biological projects involving sequencing or image data, it can be easiest to take the outputs from this pipeline as the inputs for your downstream and less computationally intensive analysis.
	Document how to run the upstream reproducible pipeline in your down steam analysis.
	This way anyone, with access to appropriate compute resources, could download your data and run the same upstream analysis to get to the same staring point for the lighter downstream analysis and all the information needed to do this is documented in the downstream analysis.
	</li>
</ul>
</details>

## 💽Environment Management / Portability

*How can people get specific versions of your software running on their systems?*

- [ ] Computational environment description provided
	- [ ]  	🥉Bronze *(Easy)*: List of package versions, e.g. output of `sessionInfo()` in R, not in a machine readable format
	- [ ]  	🥈Silver *(Intermediate)*: Structured language specific environment decription, language environment can be re-created e.g. `renv.lock` in a mostly automated fashion
	- [ ]  	🥇Gold *(Hard)*: Structured full environment description, automated ability to recreate the complete environment including system dependencies
	- [ ]  	🏆Platinum *(MAXIMUM OVERKILL)*: Your description allows the automated bootstrap of the entire* depencency tree of your environment from source with bitwise binary reproducibility (currently almost impossible to achieve, basically only approachable in Guix)

<details>
<ul>
	<li>
	The record of a specific analysis is the case where providing a complete specification of the computational environment in which code was run is perhaps the most important.
	Doing this and providing the information necessary to initiate a re-run of the analysis, in that environment, is the computational equivalent of providing protocol level methodological detail of how a bench experiment was performed.
	In addition the provision of the source of the data on which the analysis was performed and a means of both retrieving a copy of it and demonstrating that it is the same as the original input, such as hashes of the data files, is analogous to being able to get access to the same reagents and types of biological samples in an experiment.
	</li>
	<li>
	You might not need this level of detail to install a working version of a piece of software for general use.
	But for the record of a specific analysis it is ideal if we can re-run everything exactly with all the same versions if we are ever looking back to identify a potential source of error.
	</li>
	<li>
	The most approachable tools to specify and reproducibly share an interactive compute environment in which an analysis was performed are probably [binder](https://binderhub.readthedocs.io/en/latest/index.html), and [renku](https://renkulab.io/).
	</li>
	<li>
	Use of a 'lock file' (different tools may have different names for these) which specifies which software and in which versions to install in order to recreate your compute environment is ideal for this application.
	Various package and environment managers support this such as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), [renv](https://rstudio.github.io/renv/articles/renv.html), [python virtual environments](https://docs.python.org/4/library/venv.html), [poetry](https://python-poetry.org/), [nix flakes](https://zero-to-nix.com/concepts/flakes), [guix manifests](https://guix.gnu.org/en/manual/devel/en/html_node/Writing-Manifests.html).
	</li>
	<li>
	The scope of the environment managed by these tools can vary, renv for example only manage R packages, conda manages essentially any software but not necessarily all system dependencies, nix and guix can, and do, specify entire operating systems.
	</li>
	<li>
	It is common to use the more narrowly scoped package and environment management tools in conjunction with container or virtual machine build descriptions such as a Dockerfile.
	You start with an image of an operating system, add instructions to install any system dependencies and then use the environment manager to install non-system depencencies for your project into the image.
	This is incomplete as it leaves some things outside of the managed environment, it does not capture how to reconstruct the base image, and system dependencies will not necessarily be versioned.
	This is however likely the most familiar feeling experience, as it is essentially the same as what you'd do when setting up your working environment on a new computer.
	This approach covers the vast majority of cases and is a good practice to adopt now.
	Nix flakes ([vm](https://nix.dev/tutorials/nixos/nixos-configuration-on-vm.html) [container](https://nix.dev/tutorials/nixos/building-and-running-docker-images.html)) and Guix manifests ([vm](https://guix.gnu.org/manual/en/html_node/Running-Guix-in-a-VM.html) [container](https://guix.gnu.org/manual/en/html_node/Invoking-guix-pack.html)) can just generate such images directly with all system and project dependencies explicitly specified in lock files, but have a steep learning curve for their new, unfamiliar, way of working.
	They also may not yet have all the specialist software that you need packaged in their repositories.
	They can be a lot of work if they don't already have everything that you need at present but are worth watching as, once refined, can solve many of the points of friction with current approaches.
	</li>
</ul>
</details>

## 🌱 Energy Efficiency

*How can you and your users minimise wasted energy?*

- [ ] Consideration has been given to the energy efficiency of the code
	- [ ]  	🥉Bronze *(easy)*: minimise unnecessary output files
	- [ ]  	🥈Silver *(easy)*: bronze plus: Profile your code and refactor inefficient parts
	- [ ]  	🥇Gold *(intermediate)*: silver plus: Estimate and share the carbon footprint of your computations with a tools such as [green algorithms calculator](http://calculator.green-algorithms.org/)
	- [ ]  	🏆Platinum *(intermediate)*: gold plus: Offload suitable computations to hardware accelerators where possible

<details>

One off analysis code is not particularly high impact to make more efficient as it is only run a small number of times.
However, it is worth giving some consideration to the efficiency of the tools that your one off analysis might make use of or depend on.
For the most part analysis might represent a first step implementing a new method for the first time where it's correctness and comprehensibility is more important than the efficiency with which it is implemented, optimisation comes later.

Consider what can you do to make your code a little more efficient:

<ul>
	<li>
	Don't generate unnecessary outputs that will sit on people's drives unused, clean results of intermediate steps.
	You might have varying degrees of verbosity of output with a more verbose mode for debugging but defaulting to just the essentials.
	</li>
	<li>
	Good documentation and good error handing/messages can reduce the number of times people make mistakes using your code that means they re-run or partially re-run it fewer times before they figure out how to use it correctly.
	</li>
	<li>
	Some records of this type may use pipeline managers, [targets](https://books.ropensci.org/targets/) for example can integrate nicely with literate programming outputs to cache computationally expensive results, let you iterate quicky on a manuscript using those outputs and ensure that any code that needs to be re-run following a change is, all whilst being able to re-run your entire workflow and regenerate your manuscript with a single command. For pipelines in particular caching results and avoiding needing to re-compute things if possible is a good way to make best use of these features in pipeline managers for example by having small granular tasks to minimise repeated work on run failure.
	</li>
	<li>
	Choice of libraries and frameworks - some libraries may be more efficient that others or be a wrapper around an efficient implementation in another language, or be able to make use of offload to hardware accelerators.
	</li>
	<li>
	Offload to harware accelerators where available, vector matrix and array arithmathic can often benefit from very substancial speed-ups on hardware specialised for these types of calculations, or even binaries compiled with the right instruction set extensions enabled to take full advantage of hardware acceleration features on many CPUs. Doing this directly can be quite challenging but using libraries cabable of managing this offload for you can make it more approachable. (This can potentially introduce interesting reproducible computation challenges due to things like differences in handling of floating point arithmathic between hardware/firmware implementations.)
	</li>
	<li>
	Benchmarking & Profiling to locate and improve inefficient code.
	Don't optimise prematurely - it is often surpising which pieces of your code turnout to be slow, measure it first and check where to focus your attention.
	This can go hand in hand with having done robust testing as a good test suite means that you can confidently re-factor an inefficient piece of code without fear of introducing errors.
	</li>
	<li>
	Language Choice - whilst sometimes worth considering energy efficiency is rarely high up on the list of reasons to pick a programming language in this context. Familiarity of both you and anyone who might use the code after you or indeed the familiarity of the acacdemic community that might consume your code is often paramount as this maximises the speed with which you can develop your solution and others comprehend it.
	</li>
	<li>
	The people in the [Green Algorithms](https://www.green-algorithms.org/) community of practice have some useful advice and resources for anyone interersted in this subject.
	</li>
	<li>
	[The Environmentally Sustainable Computational Science forum](https://forum.escs-community.org) is a great place to find people to talk about these issues.
	</li>
</ul>
</details>


## ⚖ Governance, Conduct, & Continuity

*How can you be excellent to each other, make good decisions well, and continue to do so?*

- [ ] The project has a suitable governance model
	- [ ]  	🥉Bronze *(easy)*: The governance model is clearly communicated
	- [ ]  	🥈Silver *(easy)*: Bronze Plus - Project has continuity planning in place (2 or more from)
		- [ ] Source archived and/or mirrored to other platforms
		- [ ] Public archives of key project governance documentation and plans for continuity of operations in the events such as the loss of key project infrastructure
		- [ ] Plans of action in the event project admin(s) are no longer available
	- [ ]  	🥇Gold *(intermediate)*: Project has a governance model appropriate to its scale and goals
		- [ ] Project has clear and transparent processes
	- [ ]  	🏆Platinum *(Hard)*: Project has a track record of good governance and policy, any from:
		- [ ] Decisions have involved the appropriate person(s) and been well documented
		- [ ] Disputes are largely resolved in a respectful and amicable fashion
		- [ ] The project leadership has learned from any mistakes and implemented policy changes as a result

<details>
In the case of a record of a specific analysis it is likely that following the conclusion of that project it will become largely dormant.
You may however still encounter people refering to and making use of parts of the work and asking questions about it.
Therefore basis policies for any community spaces and user spuuort may still be sensible to indicate.
Continuity plannig is also still important as others may need to pick up where you left off at some point in the future.
<ul>
	<li>
	**Continuity planning: What happens to your project if something happens to you?**
	The code will likely live on due the distributed nature of git but what about the issue tracker, the website etc.
	Who else has the highest level of privilege on your project or a mechanism to attain it?
	The principle of least privilege dictates that you keep the number of people with this level of access to a minimum but you may then create a single point of failure.
	Password managers like [bitwarden](https://bitwarden.com/) have a feature where designated people can be given access to your vault if they request it and you do not deny it within a certain time-frame.
	This could provide a lower level admin with a mechanism to escalate their privileges if you are unable to do this for them.
	However, this delay might be an issue for continuity of operations if administrator action is needed within the waiting period.
	Game it out, have a plan, write it down, let people know you have a plan.
	</li>
	<li>
	User support
		<ul>
		<li>What support can users expect, or not expect?</li>
		<li>Where can they ask for it?</li>
		<li>Is there somewhere where users can provide support to other members of the user community, such as a forum?</li>
		<li>Can they pay for more support?</li>
		<ul>
	</li>
	<li>Whole project data longevity - what plans do you have in place to backup and archive materials pertaining to your project that are not under source control? e.g. issues in your bug tracker</li>
	<li>
	If you are the Benevolent Dictator For Life (BDFL) of your project and the Code of Conduct (CoC) is "Don't be a Dick" that's fine, for many individual hobby projects this a functional reality.
	Becoming a BDFL tends to be the default unless you take steps to avoid it and cultivate community governance as your project begins to grow - failing to do this and being stuck in charge can become quite the burden in sucessful projects.
	Be waring of adopting policies if you lack resources, time, interest, skill, or inclination to be an active enforcer, mediator and moderator of community norms and disputes, It is helpful to be clear about what you can and cannot commit to doing.
	Only by communicating this might you be able to find community members to help you with setting and enforcing these norms, if or when your community attains a scale where this becomes relevant - community management is its own skill set.
	If you can't moderate them avoid creating and/or continuing ungoverned community spaces that can become a liability for you and your project's reputation.
	Just as there are off-the-shelf licenses there are off-the-shelf codes of conduct, the [Contributor Covenant](https://www.contributor-covenant.org/) is perhaps the best known and most widely used, though may need some customisation to your needs.
	Adopting such a CoC gives you some guidance to follow if there is bad behaviour in your project's community and communicates that you as the project leadership take the responsibility of creating a respectful environment for collaboration seriously.
	It can also signal that your project is a place where everyone is wellcome but expected to treat one another with respect, and that failing to do so will result in penalties potentially including exclusion from the community.
	[The Turing Way](https://doi.org/10.5281/zenodo.3233853) provides quite a nice example of a [CoC developed specifically for their project](https://web.archive.org/web/20240412122958/https://book.the-turing-way.org/)
	You will need to provide contact information for the person(s) responsible for the enforcement of the CoC in the appropriate place and be able to follow up in the event it is used.
	git forges often recognise files with the name `CODE_OF_CONDUCT.md` in the root of project and provide a link to them on project home pages, so this is a good place to document such policies.
	If you are the BDFL of a small project then interpretation and enforcement of such a CoC tends to fall solely on you - game out some courses of action for what you'd do if faced with some common moderation challenges.
		<ul>
			<li>
				Once a project attracts a larger community there is greater scope for disputes and therefore for the need for dispute resolution mechanisms.
				Free/Libre and Open Source Software development and maintenance can be thought of as a commons so I would refer you to the [work of Elinor Ostrom](https://archive.org/details/ElinorOstromGoverningTheCommons) on how commons have been successfully (or unsuccessfully) governed when thinking about what processes to adopt for your project.
				More recently [Nathan Schneider's Governable Spaces: Democratic Design for Online Life](https://doi.org/10.1525/luminos.181) tackles some of these issues as applied to online spaces.
			</li>
			<li>
				This is summarised in the [8 Principles for Managing a Commons](https://web.archive.org/web/20240411163813/https://onthecommons.org/magazine/elinor-ostroms-8-principles-managing-commmons/index.html)
				<ol>
					<li>Define clear group boundaries.</li>
					<li>Match rules governing use of common goods to local needs and conditions.</li>
					<li>Ensure that those affected by the rules can participate in modifying the rules.</li>
					<li>Make sure the rule-making rights of community members are respected by outside authorities.</li>
					<li>Develop a system, carried out by community members, for monitoring members’ behavior.</li>
					<li>Use graduated sanctions for rule violators.</li>
					<li>Provide accessible, low-cost means for dispute resolution.</li>
					<li>Build responsibility for governing the common resource in nested tiers from the lowest level up to the entire interconnected system.</li>
				</ol>
			</li>
			<li>
				An informal [do-ocracy](https://web.archive.org/web/20230201164442/https://communityrule.info/modules/do-ocracy/) in the fiefdom of BDFL is often the default state of projects that have not given much conscious thought to how they want to be governed and are thus often subject to many of the same common failure modes of this model.
				How are decisions made in your project? Do you need the mechanisms of governance used by community and civil society organisations?
				By-laws, a committee and/or working groups, general meetings, votes, minutes? A version of these may be necessary to avoid [The Tyranny of Structurelessness](https://web.archive.org/web/20240406094310/https://www.jofreeman.com/joreen/tyranny.htm)
				How can you map these onto your development infrastructure and make the decisions of your governing bodies enactable and enforceable?
			</li>
		</ul>
	</li>
	<li>
	Does your project take donations?
	Does it have a trademark?
	Does it need a legal entity to hold these?
	Who is on the paperwork and who has signing authority?
	Who keeps track of expenditures?
	Tools & Organisations like [OpenCollective](https://opencollective.com/) can help with some of these issues.
	</li>
	<li>
	If your project has potential cybersecurity implications what procedures do you have in place for people to disclose vulnerabilities in the project so that they can be patched before they are made public.
	What systems do you have in place to disclose a vulnerability once it has been patched and ensure that users know that they need to update.
	</li>
</ul>

</details>

 <p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://gitlab.com/HDBI/data-management/checklists">Research Software Sharing, Publication, & Distribution Checklists</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://orcid.org/0000-0002-2574-9611">Richard J. Acton</a> is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>
