import React from 'react';
import { Link } from 'react-router-dom';
import './App.css';
import headerImage from './assets/chem2.png';

// Pages
import BalanceEquation from './pages/BalanceEquation';
import MolarityCalculator from './pages/MolarityCalculator';
import PHCalculator from './pages/PHCalculator';
import TitrationCurve from './pages/TitrationCurve';
import GibbsFreeEnergy from './pages/GibbsFreeEnergy';

function App() {
  const categories = [
    {
      title: "Estequiometria",
      tools: [
        { id: 'balance', name: 'Balancear Equação', path: '/balance' },
        { id: 'molarity', name: 'Calculadora de Molaridade', path: '/molarity' }
      ]
    },
    {
      title: "Equilíbrio Químico",
      tools: [
        { id: 'ph', name: 'Calculadora de pH', path: '/ph' },
        { id: 'titration', name: 'Curva de Titulação', path: '/titration' }
      ]
    },
    {
      title: "Termodinâmica",
      tools: [
        { id: 'gibbs', name: 'Energia Livre de Gibbs', path: '/gibbs' }
      ]
    }
  ];

  return (
    <div className="ferramentas-quimica">
      <div className="header-image-container">
        <img src={headerImage} alt="ChemTools Header" className="header-image" />
      </div>
      <h1>ChemTools</h1>
      <p className="description">Ferramentas essenciais para cálculos químicos</p>
      
      <nav className="container-categorias">
        <div className="button-grid">
          {categories.map(category => (
            <div key={category.title} className="category-column">
              <h2>{category.title}</h2>
              {category.tools.map(tool => (
                <Link key={tool.id} to={tool.path} className="aba">
                  {tool.name}
                </Link>
              ))}
            </div>
          ))}
        </div>
      </nav>

      <hr />
    </div>
  );
}

export default App;
