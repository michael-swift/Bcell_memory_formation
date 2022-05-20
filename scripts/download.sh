# Downloads the data from sherlock to local using rsync
# excludes file patterns I don't want
rsync -avP --exclude "*.bam*" --exclude "_*" --exclude "*SC_MULTI*" mswift2@login.sherlock.stanford.edu:/oak/stanford/groups/quake/shared/swiftShare/ilm1_processed/* .

