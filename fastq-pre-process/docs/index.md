---
layout: default
title: Home
---

This repository contains a collection of [BioWDL](https://github.com/biowdl)
workflows which can be used for aligning sequencing data. There is currently
one workflow available:
- [align-star.wdl](./align-star.html): Uses STAR to align RNA.
- [align-hisat2.wdl](./align-hisat2.html): Uses HISAT2 to align RNA.

These workflows are part of [BioWDL](https://biowdl.github.io/)
developed by the SASC team at [Leiden University Medical Center](https://www.lumc.nl/).

### Dependency requirements and tool versions
Biowdl pipelines use docker images to ensure  reproducibility. This
means that biowdl pipelines will run on any system that has docker
installed. Alternatively they can be run with singularity.

For more advanced configuration of docker or singularity please check
the [cromwell documentation on containers](
https://cromwell.readthedocs.io/en/stable/tutorials/Containers/).

Images from [biocontainers](https://biocontainers.pro) are preferred for
biowdl pipelines. The list of default images for this pipeline can be
found in the default for the `dockerImages` input.

## Contact
<p>
  <!-- Obscure e-mail address for spammers -->
For any questions about running these workflows and feature request (such as
adding additional options), please use the
<a href='https://github.com/biowdl/aligning/issues'>github issue tracker</a>
or contact the SASC team directly at: 
<a href='&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;'>
&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;</a>.
</p>
