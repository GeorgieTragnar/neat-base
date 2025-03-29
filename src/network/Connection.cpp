
#include "Connection.hpp"
#include "Node.hpp"

namespace neat {
namespace network {

int32_t Connection::getSourceId() const
{
	return source->getId();
}

int32_t Connection::getTargetId() const
{
	return target->getId();
}

double Connection::getRawInput() const
{
	return target->getValue() * weight;
}


}
}
