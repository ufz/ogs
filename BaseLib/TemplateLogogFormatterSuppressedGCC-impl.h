/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace BaseLib {

#ifdef USE_MPI
template <int T_SUPPPRESS_TOPIC_FLAG>
TemplateLogogFormatterSuppressedGCC<T_SUPPPRESS_TOPIC_FLAG>
::TemplateLogogFormatterSuppressedGCC(MPI_Comm mpi_comm)
{
    int rank = 0;
    MPI_Comm_rank(mpi_comm, &rank);
    int size = 0;
    MPI_Comm_size(mpi_comm, &size);
    int digits = size/10 + 1;
    _str_mpi_rank = "rank " + padLeft(std::to_string(rank), digits) + "| ";
}
#endif

template <int T_SUPPPRESS_TOPIC_FLAG>
LOGOG_STRING &
TemplateLogogFormatterSuppressedGCC<T_SUPPPRESS_TOPIC_FLAG>
::Format( const logog::Topic &topic, const logog::Target &target )
{
    TOPIC_FLAGS flags;
    flags = GetTopicFlags( topic );

    m_sMessageBuffer.clear();

#ifdef USE_MPI
    m_sMessageBuffer.append(_str_mpi_rank.c_str());
#endif

    if ( flags & TOPIC_FILE_NAME_FLAG )
    {
        m_sMessageBuffer.append( topic.FileName() );
        m_sMessageBuffer.append( ':' );
    }

    if ( flags & TOPIC_LINE_NUMBER_FLAG )
    {
        m_sIntBuffer.assign( topic.LineNumber() );
        m_sMessageBuffer.append( m_sIntBuffer );

        m_sMessageBuffer.append( LOGOG_CONST_STRING(": "));
    }

    RenderTimeOfDay();

    if ( flags & TOPIC_LEVEL_FLAG )
    {
        m_sMessageBuffer.append( ErrorDescription( topic.Level()));
        m_sMessageBuffer.append( LOGOG_CONST_STRING(": "));
    }

    if ( flags & TOPIC_GROUP_FLAG )
    {
        m_sMessageBuffer.append( LOGOG_CONST_STRING("{") );
        m_sMessageBuffer.append( topic.Group() );
        m_sMessageBuffer.append( LOGOG_CONST_STRING("} ") );
    }

    if ( flags & TOPIC_CATEGORY_FLAG )
    {
        m_sMessageBuffer.append( LOGOG_CONST_STRING("["));
        m_sMessageBuffer.append( topic.Category() );
        m_sMessageBuffer.append( LOGOG_CONST_STRING("] "));
    }

    if ( flags & TOPIC_MESSAGE_FLAG )
    {
        m_sMessageBuffer.append( topic.Message() );
        m_sMessageBuffer.append( (LOGOG_CHAR)'\n' );
    }

    if ( target.GetNullTerminatesStrings() )
        m_sMessageBuffer.append((LOGOG_CHAR)'\0');

    return m_sMessageBuffer;
}

} // namespace BaseLib
