import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import '../App.css';

const GibbsFreeEnergy = () => {
  const [equation, setEquation] = useState('');
  const [temperature, setTemperature] = useState('');
  const [gibbsFreeEnergy, setGibbsFreeEnergy] = useState(null);

  const handleSubmit = (e) => {
    e.preventDefault();
    // Funcionalidade para calcular a energia livre de Gibbs será implementada aqui
    setGibbsFreeEnergy('A energia livre de Gibbs calculada aparecerá aqui');
  };

  return (
    <div className="ferramentas-quimica">
      <h1>Calcular Energia Livre de Gibbs</h1>
      <p className="description">Ferramenta essencial para cálculos termodinâmicos</p>
      
      <div className="conteudo">
        <div className="button-grid">
          <div className="category-column">
            <h2>Calculadora</h2>
            <form onSubmit={handleSubmit}>
              <div>
                <label htmlFor="equation">Digite a equação química:</label>
                <input
                  type="text"
                  id="equation"
                  value={equation}
                  onChange={(e) => setEquation(e.target.value)}
                  placeholder="ex., H2 + O2 = H2O"
                  className="aba input-equation"
                />
              </div>
              <div>
                <label htmlFor="temperature">Digite a temperatura (K):</label>
                <input
                  type="number"
                  id="temperature"
                  value={temperature}
                  onChange={(e) => setTemperature(e.target.value)}
                  placeholder="ex., 298"
                  className="aba input-equation"
                />
              </div>
              <button type="submit" className="aba">Calcular Energia Livre de Gibbs</button>
            </form>
            {gibbsFreeEnergy && (
              <div className="result">
                <h3>Energia Livre de Gibbs:</h3>
                <p>{gibbsFreeEnergy}</p>
              </div>
            )}
          </div>
          <div className="category-column">
            <h2>Descrição do Problema</h2>
            <p>
              A energia livre de Gibbs é uma medida termodinâmica fundamental que determina
              a espontaneidade de uma reação química. Esta ferramenta ajuda você a calcular
              a energia livre de Gibbs para uma dada reação química em uma temperatura específica.
            </p>
            <p>
              Para usar, insira a equação química no primeiro campo e a temperatura em Kelvin
              no segundo campo. Clique em "Calcular Energia Livre de Gibbs" e a ferramenta
              irá determinar se a reação é espontânea (ΔG {'<'} 0), em equilíbrio (ΔG = 0),
              ou não espontânea (ΔG {'>'} 0) nas condições especificadas.
            </p>
          </div>
        </div>
      </div>
      <Link to="/" className="aba">Voltar para a Página Inicial</Link>
    </div>
  );
};

export default GibbsFreeEnergy;
