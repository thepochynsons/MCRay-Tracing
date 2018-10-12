#include "inputmanager.h"
#include <SFML/Graphics.hpp>

inputManager::inputManager() :
    translationSpeedX(0.001f),
    translationSpeedY(0.001f),
    translationSpeedZ(0.001f)
{

}


void inputManager::setTransducer(transducer_ * transducer)
{
    this->transducer = transducer;
}

void inputManager::update(float elapsedTimeInMillis)
{
    assert(transducer != nullptr);
    for (sf::RenderWindow * w : windows)
    {
        if (!w->isOpen())
        {
            //DEBUGMSG("Window not open. Input handling could not be done.");
            std::cout<< "error";
            return;
        }

        sf::Event event;
        while (w->pollEvent(event)) {
            if (event.type == sf::Event::Closed)
            {
                isQuitPressed = true;
                exit(0);
            }
            else if (event.type == sf::Event::Resized)
            {
                w->setView(sf::View(sf::FloatRect(0, 0,
                                                 event.size.width,
                                                 event.size.height)));
            }
        }
    }

      if (inputEnabled)
    {
        updateKeyboard(elapsedTimeInMillis);

        transducer->update();
    }
}



void inputManager::updateKeyboard(float elapsedTimeInMillis)
{
    handleKeyboardMovement(elapsedTimeInMillis);
}

void inputManager::handleKeyboardMovement(float elapsedTimeInMillis)
{
    float rx = 0.0f, ry = 0.0f, rz = 0.0f;
    float dx = 0.0f, dy = 0.0f, dz = 0.0f;
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::R))
    {
   //     transducer->resetPosition();// si hacemos esta funcion seria que vuelva a partir de la posicion del archivo de configuracion
    }
    bool translationDetected = false;

    // translate Y
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::Add) ||
        sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
    {
        //dy = 100 * translationSpeedY * elapsedTimeInMillis;
        dy = 0.1;
        translationDetected = true;
    }
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::Subtract) ||
        sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
    {
        //dy = -100 * translationSpeedY * elapsedTimeInMillis;
        dy = -0.1;
        translationDetected = true;
    }

    // translate Z
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::Num1))
    {
        //dz = 100 * translationSpeedZ * elapsedTimeInMillis;
        dz = 0.1;
        translationDetected = true;
    }
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::Num2))
    {
        //dz = -100 * translationSpeedZ * elapsedTimeInMillis;
        dz = -0.1;
        translationDetected = true;
    }

    // translate X
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad4) ||
        sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
    {
        //dx = 100 * translationSpeedX * elapsedTimeInMillis;
        dx = 0.1;
        translationDetected = true;
    }
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad6) ||
        sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
    {
        //dx = -100 * translationSpeedX * elapsedTimeInMillis;
        dx = -0.1;
        translationDetected = true;
    }

    if (translationDetected)
    {
        btVector3 lastPosition = transducer->getPosition();
        transducer->setPosition(btVector3(lastPosition[0] + dx, lastPosition[1] + dy, lastPosition[2] + dz));
    }
}



