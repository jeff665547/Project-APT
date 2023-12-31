////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////
#include "rtest/RT_Args.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Args::RT_Args()
{
  setName("");
  setValue("");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Args::RT_Args(const std::string &name, const std::string &value)
{
  setName(name);
  setValue(value);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Args::~RT_Args()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Args::dump()
{
  Verbose::out(2, "***RT_Args Dump***");
  Verbose::out(2, "Name: " + this->m_Name);
  Verbose::out(2, "Value: " + this->m_Value);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Args::setName(const std::string &name)
{
  this->m_Name = name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Args::setValue(const std::string &value)
{
  this->m_Value = value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string RT_Args::getName()
{
  return this->m_Name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string RT_Args::getValue()
{
  return this->m_Value;
}
