
DIRS = evalresp libmseed src

# Test for Makefile/makefile and run make, run configure if it exists
# and no Makefile does.

# As a special case for evalresp do not pass targets except "clean".

all clean static install gcc gcc32 gcc64 debug gccdebug gcc32debug gcc64debug ::
	@for d in $(DIRS) ; do \
	  if [ ! -f $$d/Makefile -a ! -f $$d/makefile ] ; then \
	    if [ -x $$d/configure ] ; then \
	      echo "Running configure in $$d" ; \
	      ( cd $$d && ./configure --enable-lib-mode --enable-log-label ) ; \
	    fi ; \
	  fi ; \
	  echo "Running $(MAKE) $@ in $$d" ; \
	  if [ -f $$d/Makefile -o -f $$d/makefile ] ; then \
	    if [ "$$d" = "evalresp" -a "$@" != "clean" ] ; then \
	      ( cd $$d && $(MAKE) ) ; \
	    else \
	      ( cd $$d && $(MAKE) $@ ) ; \
	    fi ; \
	  elif [ -d $$d ] ; \
	    then ( echo "ERROR: no Makefile/makefile in $$d for $(CC)" ) ; \
	  fi ; \
	done

