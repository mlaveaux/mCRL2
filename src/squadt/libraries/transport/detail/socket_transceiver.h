#ifndef SOCKET_TRANSCEIVER_H
#define SOCKET_TRANSCEIVER_H

#include <boost/bind.hpp>
#include <boost/asio.hpp>

#include <transport/transporter.h>
#include <transport/detail/transceiver.h>
#include <transport/detail/socket_scheduler.h>

namespace transport {
  namespace transceiver {

    /* Class that is used internally for direct transmitting/receiving */
    class socket_transceiver : public basic_transceiver {
      friend class transport::transporter;
      friend class transport::listener::socket_listener;
  
      public:
        /** \brief IP version 4 address verifier (refer to the asio documentation) */
        typedef asio::ipv4::address address;

        /** \brief IP version 4 host class (refer to the asio documentation) */
        typedef asio::ipv4::host    host;

        /** \brief Convenience type to hide the boost shared pointer implementation */
        typedef boost::shared_ptr < socket_transceiver > ptr;

      private:

        /** \brief Host name resolver */
        static asio::ipv4::host_resolver resolver;

        /** \brief Scheduler for asynchronous socket communication */
        static socket_scheduler          scheduler;

        /** \brief Default port for socket connections */
        static long                      default_port;

        /** \brief Size of the input buffer */
        static unsigned int              input_buffer_size;

        /** \brief The input buffer */
        char*                            buffer;

        /** \brief The local endpoint of a connection */
        asio::stream_socket              socket;
 
      private:

        /** \brief Send a string input stream to the peer */
        void send(const std::string& data);
  
        /** \brief Send the contents of an input stream to the peer */
        void send(std::istream& data);

        /** \brief Connect to a peer using an address and a port */
        void connect(const address& address = address::loopback(), const long port = default_port);

        /** \brief Connect to a peer using an address and a port */
        void connect(const std::string& host_name, const long port = default_port);

        /** \brief Returns an object with the local hosts name and addresses */
        static host get_local_host();

        /** \brief Read from the socket */
        void handle_receive(const asio::error& e);

        /** \brief Process results from a write operation on the socket */
        void handle_write(const asio::error& e);

        /** \brief Start listening for new data */
        void activate();

      public:

        /** \brief Constructor that connects to a port on an address */
        socket_transceiver(transporter& o);

        /** \brief Terminate the connection with the peer */
        void disconnect(transporter::connection_ptr);

        /** \brief Destructor */
        ~socket_transceiver();
    };

    inline void socket_transceiver::disconnect(transporter::connection_ptr) {
      socket.shutdown(asio::stream_socket::shutdown_both);
      socket.close();

      basic_transceiver::handle_disconnect(this);
    }
  }
}

#endif

