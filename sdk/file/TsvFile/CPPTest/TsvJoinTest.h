////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


/// @file   TsvJoinTest.h
/// @brief  Header for TsvJoinTest.cpp

#ifndef TSV_JOIN_TEST_H
#define TSV_JOIN_TEST_H

#include <cppunit/extensions/HelperMacros.h>

class TsvJoinTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TsvJoinTest );
  CPPUNIT_TEST( testTsvJoin );
  CPPUNIT_TEST_SUITE_END();

public:
  void testTsvJoin();
};

#endif /* TSV_JOIN_TEST_H */
