#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"
#include "display_tools.h"

int main(){

    vector<Matrix> localKs;
    vector<Vector> localbs;

    Matrix K;
    Vector b;
    Vector T;

    mesh m;
    leerMallayCondiciones(m);

    crearSistemasLocales(m,localKs,localbs);

    zeroes(K,m.getSize(NODES));
    zeroes(b,m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);

    applyNeumann(m,b);

    applyDirichlet(m,K,b);
    zeroes(T,b.size());
    calculate(K,b,T);

    cout << "La respuesta es: \n";
    showVector(T);

    return 0;
}
