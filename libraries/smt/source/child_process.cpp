// Author(s): Thomas Neele
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file child_process.cpp

#include "mcrl2/smt/child_process.h"
#include "mcrl2/utilities/platform.h"

#ifdef MCRL2_PLATFORM_WINDOWS
  #include <windows.h>
  #include <tchar.h>
  #include <cstdio>
  #include <strsafe.h>
#else
  #include <unistd.h>
  #include <sys/select.h>
  #include <sys/wait.h>
#endif // MCRL2_PLATFORM_WINDOWS

#include <cstring>
#include <csignal>
#include <cerrno>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::smt
{

#ifdef MCRL2_PLATFORM_WINDOWS

struct child_process::platform_impl
{
  HANDLE g_hChildStd_IN_Rd = NULL;
  HANDLE g_hChildStd_IN_Wr = NULL;
  HANDLE g_hChildStd_OUT_Rd = NULL;
  HANDLE g_hChildStd_OUT_Wr = NULL;
};

void child_process::write(const std::string& command) const
{
  DWORD dwWritten;
  BOOL bSuccess = FALSE;

  bSuccess = WriteFile(m_pimpl->g_hChildStd_IN_Wr, command.c_str(), static_cast<DWORD>(command.size()), &dwWritten, NULL);
  if(!bSuccess || dwWritten != command.size())
  {
    throw mcrl2::runtime_error("Failed to write SMT problem to the solver.\n" + command);
  }
}


std::string child_process::read() const
{
  DWORD dwRead, totalRead;
  CHAR output[512];
  BOOL bSuccess = FALSE;

  totalRead = 0;
  do
  {
    bSuccess = ReadFile( m_pimpl->g_hChildStd_OUT_Rd, output + totalRead, 512 - totalRead, &dwRead, NULL);
    totalRead += dwRead;
  } while(bSuccess && output[totalRead-1] != '\n');
  if(!bSuccess || totalRead == 0 )
  {
    throw mcrl2::runtime_error("Reading from child process " + m_name + " failed.");
  }
  return std::string(output, totalRead);
}

std::string child_process::read(const std::chrono::microseconds& timeout) const
{
  // Not implemented on Windows
  return read();
}

void child_process::send_sigint() const
{
  // Not implemented on Windows
}

void child_process::initialize()
{
  m_pimpl = std::make_shared<child_process::platform_impl>();
  SECURITY_ATTRIBUTES saAttr;
  // Set the bInheritHandle flag so pipe handles are inherited.
  saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
  saAttr.bInheritHandle = TRUE;
  saAttr.lpSecurityDescriptor = NULL;

  // Create a pipe for the child process's STDOUT.
  if (!CreatePipe(&m_pimpl->g_hChildStd_OUT_Rd, &m_pimpl->g_hChildStd_OUT_Wr, &saAttr, 0))
  {
    throw mcrl2::runtime_error("Could not connect to SMT solver: failed to create pipe Stdout");
  }
  // Ensure the read handle to the pipe for STDOUT is not inherited.
  if(!SetHandleInformation(m_pimpl->g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0))
  {
    throw mcrl2::runtime_error("Could not modify SMT solver handle: SetHandleInformation Stdout");
  }
  // Create a pipe for the child process's STDIN.
  if(!CreatePipe(&m_pimpl->g_hChildStd_IN_Rd, &m_pimpl->g_hChildStd_IN_Wr, &saAttr, 0))
  {
    throw mcrl2::runtime_error("Could not connect to SMT solver: failed to create pipe Stdin");
  }
  // Ensure the write handle to the pipe for STDIN is not inherited.
  if(!SetHandleInformation(m_pimpl->g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0))
  {
    throw mcrl2::runtime_error("Could not modify SMT solver handle: SetHandleInformation Stdin");
  }
  // Create the child process.
  TCHAR szCmdline[] = TEXT("z3.exe -smt2 -in");
  PROCESS_INFORMATION piProcInfo;
  STARTUPINFO siStartInfo;
  BOOL bSuccess = FALSE;

  // Set up members of the PROCESS_INFORMATION structure.
  ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));

  // Set up members of the STARTUPINFO structure.
  // This structure specifies the STDIN and STDOUT handles for redirection.
  ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
  siStartInfo.cb = sizeof(STARTUPINFO);
  siStartInfo.hStdError = m_pimpl->g_hChildStd_OUT_Wr;
  siStartInfo.hStdOutput = m_pimpl->g_hChildStd_OUT_Wr;
  siStartInfo.hStdInput = m_pimpl->g_hChildStd_IN_Rd;
  siStartInfo.dwFlags |= STARTF_USESTDHANDLES;

  // Create the child process.
  bSuccess = CreateProcess(NULL,
    szCmdline,     // command line
    NULL,          // process security attributes
    NULL,          // primary thread security attributes
    TRUE,          // handles are inherited
    0,             // creation flags
    NULL,          // use parent's environment
    NULL,          // use parent's current directory
    &siStartInfo,  // STARTUPINFO pointer
    &piProcInfo);  // receives PROCESS_INFORMATION

  // If an error occurs, exit the application.
  if(!bSuccess)
  {
    throw mcrl2::runtime_error("Problem occurred while creating SMT solver process.\nIs " + m_name + " in your path?");
  }
  else
  {
    // Close handles to the child process and its primary thread.
    // Some applications might keep these handles to monitor the status
    // of the child process, for example.
    CloseHandle(piProcInfo.hProcess);
    CloseHandle(piProcInfo.hThread);
  }
}

child_process::~child_process()
{
  // Clean up the child process
  // It terminates when we close its stdin
  CloseHandle(m_pimpl->g_hChildStd_IN_Wr);
}


#else // MCRL2_PLATFORM_WINDOWS

struct child_process::platform_impl
{
  pid_t child_pid;
  int pipe_stdin[2];
  int pipe_stdout[2];
  int pipe_stderr[2];
};

void child_process::write(const std::string& command) const
{
  if(::write(m_pimpl->pipe_stdin[1], command.c_str(), command.size()) < 0)
  {
    throw mcrl2::runtime_error("Failed to write to " + m_name + ".\n" + command);
  }
}

std::string child_process::read() const
{
  char output[512];
  // check return value
  int count = ::read(m_pimpl->pipe_stdout[0], output, 512);
  if(count > 0)
  {
    return std::string(output, count);
  }
  else
  {
    std::string message;

    while(int i = ::read(m_pimpl->pipe_stderr[0], output, 512))
    {
      message.append(output, 0, i);
    }

    throw mcrl2::runtime_error("Reading from child process " + m_name + " failed: " + message);
  }
}

std::string child_process::read(const std::chrono::microseconds& timeout) const
{
  struct timeval tv{};
  tv.tv_sec = 0;
  tv.tv_usec = timeout.count();

  fd_set readfds;
  FD_ZERO(&readfds);
  FD_SET(m_pimpl->pipe_stdout[0], &readfds);

  int result = ::select(m_pimpl->pipe_stdout[0] + 1, &readfds, nullptr, nullptr, &tv);

  if(result == -1)
  {
    throw mcrl2::runtime_error("Error while waiting for input from " + m_name);
  }
  else if(result == 0)
  {
    send_sigint();
  }
  return read();
}

void child_process::send_sigint() const
{
  ::kill(m_pimpl->child_pid, SIGINT);
}

void child_process::initialize()
{
  signal(SIGPIPE, SIG_IGN);
  m_pimpl = std::make_shared<child_process::platform_impl>();
  // Create pipes (two pairs r/w)
  if (::pipe(&m_pimpl->pipe_stdin[0]) < 0)
  {
    throw mcrl2::runtime_error("failed to create pipe");
  }

  if (::pipe(&m_pimpl->pipe_stdout[0]) < 0)
  {
    throw mcrl2::runtime_error("failed to create pipe");
  }

  if (::pipe(&m_pimpl->pipe_stderr[0]) < 0)
  {
    throw mcrl2::runtime_error("failed to create pipe");
  }

  // fork process
  pid_t pid = ::fork();
  if (pid == 0)
  {
    // Forked process keeps waiting for input
    ::dup2(m_pimpl->pipe_stdin[0], STDIN_FILENO);
    ::dup2(m_pimpl->pipe_stdout[1], STDOUT_FILENO);
    ::dup2(m_pimpl->pipe_stderr[1], STDERR_FILENO);

    ::close(m_pimpl->pipe_stdin[1]);
    ::close(m_pimpl->pipe_stdout[0]);
    ::close(m_pimpl->pipe_stderr[0]);

    ::execlp("z3", "z3", "-smt2", "-in", (char*) nullptr);

    ::_exit(errno);
  }
  else if (pid < 0)
  {
    ::close(m_pimpl->pipe_stdin[0]);
    ::close(m_pimpl->pipe_stdin[1]);
    ::close(m_pimpl->pipe_stdout[0]);
    ::close(m_pimpl->pipe_stdout[1]);
    ::close(m_pimpl->pipe_stderr[0]);
    ::close(m_pimpl->pipe_stderr[1]);

    throw mcrl2::runtime_error("Problem occurred while forking SMT solver process: " + std::string(strerror(errno)) + ".\nIs " + m_name + " in your path?");
  }
  else
  {
    m_pimpl->child_pid = pid;
    ::close(m_pimpl->pipe_stdin[0]);
    ::close(m_pimpl->pipe_stdout[1]);
    ::close(m_pimpl->pipe_stderr[1]);
  }
}

child_process::~child_process()
{
  // Clean up the child process
  // It terminates when we close its stdin
  ::close(m_pimpl->pipe_stdin[1]);
  ::close(m_pimpl->pipe_stdout[0]);
  ::close(m_pimpl->pipe_stderr[0]);

  int return_status;
  ::wait(&return_status);
}

#endif // MCRL2_PLATFORM_WINDOWS


} // namespace mcrl2::smt

