/**
 * \file socket.hpp Cross-platform socket abstractions
 */

/**
 * \todo Socket support is unfinished and is not working yet.  Please review socket.hpp if you'd like to implement
 * a cross-platform socket abstraction.
 */
#ifndef __LOGOG_SOCKET_HPP_
#define __LOGOG_SOCKET_HPP_

#ifdef LOGOG_FLAVOR_WINDOWS
#pragma comment(lib, "wsock32.lib")
#endif

#include <memory.h>

namespace logog
{
/** A sending or receiving socket to be used as a target. */
class Socket : public Target
{
public:
    static int Initialize()
    {
        static bool bInitialized = false;

        if ( !bInitialized )
        {
#ifdef LOGOG_FLAVOR_WINDOWS
            WORD wVersionRequested;
            WSADATA wsaData;
            int err;
            wVersionRequested = MAKEWORD( 2, 2 );

            err = WSAStartup( wVersionRequested, &wsaData );

            if ( err != 0 )
            {
#ifdef LOGOG_INTERNAL_DEBUGGING
                LOGOG_COUT << _LG("WSAStartup failed with error: ") <<  err << endl;
#endif
                return 1;
            }

#endif
        }

        bInitialized = true;

        return 0;
    }

    static void Shutdown()
    {
#ifdef LOGOG_FLAVOR_WINDOWS
        WSACleanup();
#endif
    }

    static const int MAXHOSTNAME = 255;

    Socket(
        int type = SOCK_STREAM,
        int port = LOGOG_DEFAULT_PORT
    )
    {
        m_Socket = -1;
        m_nType = type;
        m_nPort = port;
    }

    virtual void Close()
    {
#ifdef LOGOG_FLAVOR_WINDOWS
        closesocket( m_Socket );
#endif
#ifdef LOGOG_FLAVOR_POSIX
        close( m_Socket );
#endif
    }

#ifdef NYI
    virtual int Create( int type,
                        int port )
    {
        char myname[MAXHOSTNAME+1];
        int err;
        struct sockaddr_in sa;
        struct hostent *hp;

        memset( &sa, 0, sizeof( struct sockaddr_in ) ); /* clear our address */
        gethostname( myname, MAXHOSTNAME );         /* who are we? */
        hp= gethostbyname( myname );                /* get our address info */

        if ( hp == NULL )                           /* we don't exist !? */
            return -1;

        sa.sin_family= hp->h_addrtype;              /* this is our host address */
        sa.sin_port= (u_short)htons( (u_short)port );                /* this is our port number */

        if (( m_Socket = socket( AF_INET, type, 0 ) ) < 0 ) /* create socket */
            return -1;

        if ( bind( m_Socket,  (struct sockaddr*) &sa, sizeof( struct sockaddr_in ) ) < 0 )
        {
            Close();
            return -1;
        }

        if (( err = SetNonBlocking() ) != 0 )
            return err;
    }

#endif // NYI

    virtual int SetNonBlocking()
    {
        int err;

#ifdef LOGOG_FLAVOR_POSIX
        int flags;
        flags = socket_fcntl( m_Socket, F_GETFL, 0 );
        flags |= O_NONBLOCK;
        err = socket_fcntl( m_Socket, F_SETFL, flags );
#endif

#ifdef LOGOG_FLAVOR_WINDOWS
        unsigned long parg;
        parg = 1;
        err = ioctlsocket( m_Socket, FIONBIO, &parg );
#endif
        return err;
    }

    virtual int Output( const LOGOG_STRING &output ) = 0;

protected:
    int m_Socket;
    int m_nType;
    int m_nPort;
};

class SocketServer : Socket
{

};
}

#endif // __LOGOG_CHECKPOINT_HPP_
