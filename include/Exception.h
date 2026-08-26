// Exception.h
// (C) Coire Cadeau 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <exception>

class Exception : public std::exception {
  const std::string message;

 public:
  Exception(const char* msg) : message(msg) {}
  virtual ~Exception() throw() {}
  const char* what() const throw() { return message.c_str(); }
};

#endif // EXCEPTION_H
