#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2011-2014 Clemson University.
# 
#  This file was originally written by Bradley S. Meyer.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#/////////////////////////////////////////////////////////////////////////////*/

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate screen examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef NUCNET_TARGET
NUCNET_TARGET = ../..
endif

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

CC_SCREEN = $(CC) -I $(NUCNET_TARGET)

#===============================================================================
# Objects.
#===============================================================================

SCREEN_OBJ = $(WN_OBJ)        \
             $(USER_OBJ)      \
             $(NNT_OBJ) 

#===============================================================================
# Executables.
#===============================================================================

SCREEN_EXEC = compare_screened_rates_in_zones			\
              compute_screening_factors_in_zones		\
              compute_reverse_ratio_correction_in_zones		\

$(SCREEN_EXEC): $(SCREEN_OBJ)
	$(CC_SCREEN) $(SCREEN_OBJ) -o $(BINDIR)/$@ $@.cpp $(CLIBS)

.PHONY all_screen: $(SCREEN_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_screen cleanall_screen

clean_screen:
	rm -f $(SCREEN_OBJ)

cleanall_screen: clean_screen
	rm -f $(BINDIR)/$(SCREEN_EXEC) $(BINDIR)/$(THERMO_EXEC).exe
