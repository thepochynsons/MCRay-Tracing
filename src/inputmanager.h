#ifndef INPUTMANAGER_H
#define INPUTMANAGER_H

#include <vector>
#include "transducer.h"

namespace sf
{
    class RenderWindow;
}

class inputManager
{
    using transducer_ = transducer<512>;

public:
        inputManager();
        //~inputmanager();
        void setTransducer(transducer_ * transducer);
        void update(float elapsedTimeInMillis);
        void updateKeyboard(float elapsedTimeInMillis);
        void handleKeyboardMovement(float elapsedTimeInMillis);

private:
        transducer_ * transducer;
        float translationSpeedX, translationSpeedY, translationSpeedZ;
        bool inputEnabled = true, isQuitPressed = false;
        std::vector<sf::RenderWindow*> windows;


};
#endif // INPUTMANAGER_H
